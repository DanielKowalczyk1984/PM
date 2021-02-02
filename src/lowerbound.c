#include "gurobi_c.h"
#include "job.h"
#include "lp.h"
#include "math.h"
#include "pricingstabilizationwrapper.h"
#include "scheduleset.h"
#include "solver.h"
#include "util.h"
#include "wct.h"
#include "wctparms.h"
#include "wctprivate.h"

static const double min_nb_del_row_ratio = 0.9;

void g_print_ages_col(gpointer data, MAYBE_UNUSED gpointer user_data) {
    ScheduleSet* x = (ScheduleSet*)data;

    printf(" %4d", x->age);
}

/** Help function for column generation */
static void print_ages(NodeData* pd) {
    printf("AGES:");

    g_ptr_array_foreach(pd->localColPool, g_print_ages_col, NULL);

    printf("\n");
}

void g_grow_ages(gpointer data, gpointer user_data) {
    ScheduleSet* x = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;

    if (pd->column_status[x->id] == lp_interface_LOWER ||
        pd->column_status[x->id] == lp_interface_FREE) {
        x->age++;

        if (x->age > pd->retirementage) {
            pd->zero_count++;
        }
    } else {
        x->age = 0;
    }
}

static int grow_ages(NodeData* pd) {
    int val = 0;
    int nb_cols = 0;
    lp_interface_get_nb_cols(pd->RMP, &nb_cols);
    assert(nb_cols - pd->id_pseudo_schedules == pd->localColPool->len);
    CC_IFFREE(pd->column_status, int);
    if (pd->localColPool->len > 0) {
        pd->column_status = CC_SAFE_MALLOC(pd->localColPool->len, int);
        CCcheck_NULL_2(pd->column_status, "Failed to allocate column_status");
        val = lp_interface_basis_cols(pd->RMP, pd->column_status,
                                      pd->id_pseudo_schedules);
        CCcheck_val_2(val, "Failed in lp_interface_basis_cols");
        pd->zero_count = 0;

        g_ptr_array_foreach(pd->localColPool, g_grow_ages, pd);
    }

CLEAN:
    return val;
}

int delete_unused_rows(NodeData* pd) {
    int     val = 0;
    int     nb_rows = 0;
    double* slack = &g_array_index(pd->slack, double, 0);
    GArray* del_indices = g_array_new(FALSE, FALSE, sizeof(int));

    lp_interface_get_nb_rows(pd->RMP, &nb_rows);
    // assert(nb_rows == pd->slack->len);
    lp_interface_slack(pd->RMP, slack);

    int it = pd->id_valid_cuts;
    int first_del = -1;
    int last_del = -1;
    for (int i = pd->id_valid_cuts; i < nb_rows; i++) {
        if (fabs(slack[i]) < EPS) {
            if (first_del != -1) {
                val = delete_unused_rows_range(pd, first_del, last_del);
                CCcheck_val(val, "Failed in deleterows lp_interface");
                it = it - (last_del - first_del);
                first_del = last_del = -1;
            } else {
                it++;
            }
        } else {
            if (first_del == -1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            it++;
        }
    }

    if (first_del != -1) {
        delete_unused_rows_range(pd, first_del, last_del);
    }

    call_update_rows_coeff(pd);

    return val;
}

int delete_old_schedules(NodeData* pd) {
    int          val = 0;
    int          min_numdel = floor(pd->nb_jobs * min_nb_del_row_ratio);
    int          nb_col = 0;
    guint        i = 0U;
    guint        count = pd->localColPool->len;
    ScheduleSet* tmp_schedule = (ScheduleSet*)NULL;
    /** pd->zero_count can be deprecated! */
    pd->zero_count = 0;

    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        if (tmp_schedule->age > 0) {
            pd->zero_count++;
        }
    }

    if (pd->zero_count > min_numdel) {
        int it = 0;
        int first_del = -1;
        int last_del = -1;
        lp_interface_get_nb_cols(pd->RMP, &nb_col);
        assert(nb_col - pd->id_pseudo_schedules == count);
        for (i = 0; i < count; ++i) {
            tmp_schedule =
                (ScheduleSet*)g_ptr_array_index(pd->localColPool, it);
            if (tmp_schedule->age <= pd->retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = lp_interface_deletecols(
                        pd->RMP, first_del + pd->id_pseudo_schedules,
                        last_del + pd->id_pseudo_schedules);
                    CCcheck_val_2(val, "Failed in lp_interface_deletecols");
                    g_ptr_array_remove_range(pd->localColPool, first_del,
                                             last_del - first_del + 1);
                    it = it - (last_del - first_del);
                    first_del = last_del = -1;
                } else {
                    it++;
                }
            } else {
                if (first_del == -1) {
                    first_del = it;
                    last_del = first_del;
                } else {
                    last_del++;
                }
                it++;
            }
        }

        if (first_del != -1) {
            lp_interface_deletecols(pd->RMP,
                                    first_del + pd->id_pseudo_schedules,
                                    last_del + pd->id_pseudo_schedules);
            CCcheck_val_2(val, "Failed in lp_interface_deletecols");
            g_ptr_array_remove_range(pd->localColPool, first_del,
                                     last_del - first_del + 1);
        }

        if (dbg_lvl() > 1) {
            printf("Deleted %d out of %d columns with age > %d.\n",
                   pd->zero_count, count, pd->retirementage);
        }

        lp_interface_get_nb_cols(pd->RMP, &nb_col);
        assert(pd->localColPool->len == nb_col - pd->id_pseudo_schedules);
        for (i = 0; i < pd->localColPool->len; ++i) {
            tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
            tmp_schedule->id = i;
        }
        pd->zero_count = 0;
    }

CLEAN:
    return val;
}

int delete_infeasible_schedules(NodeData* pd) {
    int          val = 0;
    int          nb_col = 0;
    guint        i = 0;
    guint        count = pd->localColPool->len;
    ScheduleSet* tmp_schedule = (ScheduleSet*)NULL;
    /** pd->zero_count can be deprecated! */
    pd->zero_count = 0;

    int it = 0;
    int first_del = -1;
    int last_del = -1;
    lp_interface_get_nb_cols(pd->RMP, &nb_col);
    assert(nb_col - pd->id_pseudo_schedules == count);
    for (i = 0; i < count; ++i) {
        // while (it < pd->localColPool->len) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, it);
        if (tmp_schedule->del != 0) {
            if (first_del != -1) {
                /** Delete recently found deletion range.*/
                val = lp_interface_deletecols(
                    pd->RMP, first_del + pd->id_pseudo_schedules,
                    last_del + pd->id_pseudo_schedules);
                CCcheck_val_2(val, "Failed in lp_interface_deletecols");
                g_ptr_array_remove_range(pd->localColPool, first_del,
                                         last_del - first_del + 1);
                pd->zero_count += last_del - first_del + 1;
                it = it - (last_del - first_del);
                first_del = last_del = -1;
            } else {
                it++;
            }
        } else {
            if (first_del == -1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            it++;
        }
    }

    if (first_del != -1) {
        lp_interface_deletecols(pd->RMP, first_del + pd->id_pseudo_schedules,
                                last_del + pd->id_pseudo_schedules);
        CCcheck_val_2(val, "Failed in lp_interface_deletecols");
        g_ptr_array_remove_range(pd->localColPool, first_del,
                                 last_del - first_del + 1);
        pd->zero_count += last_del - first_del + 1;
    }

    if (dbg_lvl() > 1) {
        printf(
            "Deleted %d out of %d columns(infeasible columns after reduce cost "
            "fixing).\n",
            pd->zero_count, count);
    }

    lp_interface_get_nb_cols(pd->RMP, &nb_col);
    assert(pd->localColPool->len == nb_col - pd->id_pseudo_schedules);
    if (dbg_lvl() > 1) {
        printf("number of cols = %d\n", nb_col - pd->id_pseudo_schedules);
    }

    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        tmp_schedule->id = i;
    }

    if (pd->zero_count > 0) {
        solve_relaxation(pd);
        pd->update = 1;
    }

CLEAN:
    return val;
}

void g_make_pi_feasible(gpointer data, gpointer user_data) {
    ScheduleSet* x = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;
    Job*         tmp_j = (Job*)NULL;

    double colsum = .0;

    for (guint i = 0; i < x->job_list->len; ++i) {
        tmp_j = (Job*)g_ptr_array_index(x->job_list, i);
        if (signbit((double)pd->pi->data[tmp_j->job])) {
            pd->pi->data[tmp_j->job] = 0.0;
        }

        colsum += pd->pi->data[tmp_j->job];
        colsum = nextafter(colsum, DBL_MAX);
    }

    if (!signbit((double)pd->pi->data[pd->nb_jobs])) {
        g_array_index(pd->pi, double, pd->nb_jobs) = 0.0;
    }

    colsum += (double)pd->pi->data[pd->nb_jobs];
    colsum = nextafter(colsum, DBL_MAX);

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (guint i = 0; i < x->job_list->len; ++i) {
            tmp_j = (Job*)g_ptr_array_index(x->job_list, i);
            g_array_index(pd->pi, double, tmp_j->job) /= colsum;
            g_array_index(pd->pi, double, tmp_j->job) *=
                x->total_weighted_completion_time;
            newcolsum += g_array_index(pd->pi, double, tmp_j->job);
        }

        g_array_index(pd->pi, double, pd->nb_jobs) /= colsum;
        g_array_index(pd->pi, double, pd->nb_jobs) *=
            x->total_weighted_completion_time;
        newcolsum += g_array_index(pd->pi, double, pd->nb_jobs);

        if (dbg_lvl() > 1) {
            printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                   x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible(NodeData* pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible, pd);
}

void g_make_pi_feasible_farkas(gpointer data, gpointer user_data) {
    ScheduleSet* x = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;

    double colsum = .0;

    for (guint i = 0; i < x->job_list->len; ++i) {
        double* tmp = &g_array_index(pd->pi, double, i);
        if (signbit(*tmp)) {
            *tmp = 0.0;
        }

        colsum += *tmp;
        colsum = nextafter(colsum, DBL_MAX);
    }

    colsum += g_array_index(pd->pi, double, pd->nb_jobs);

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (guint i = 0; i < x->job_list->len; ++i) {
            double* tmp = &g_array_index(pd->pi, double, i);
            *tmp /= colsum;
            *tmp *= x->total_weighted_completion_time;
            newcolsum += *tmp;
        }

        double* tmp = &g_array_index(pd->pi, double, pd->nb_jobs);
        *tmp /= colsum;
        *tmp *= x->total_weighted_completion_time;
        newcolsum += *tmp;

        if (dbg_lvl() > 1) {
            printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                   x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible_farkas_pricing(NodeData* pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible_farkas, pd);
}

int compute_objective(NodeData* pd) {
    int val = 0;
    // int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    double* tmp = &g_array_index(pd->pi, double, 0);
    double* tmp_rhs = &g_array_index(pd->rhs, double, 0);

    for (int i = 0; i < pd->nb_rows; i++) {
        // if (i != pd->nb_jobs) {
        // pd->eta_out += tmp[i] * tmp_rhs[i];
        // }
        pd->LP_lower_bound_dual += tmp[i] * tmp_rhs[i];
    }
    pd->LP_lower_bound_dual -= EPS_BOUND;

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = lp_interface_objval(pd->RMP, &(pd->LP_lower_bound));
    pd->LP_lower_bound -= EPS_BOUND;
    CCcheck_val_2(val, "lp_interface_objval failed");
    pd->lower_bound =
        ((int)ceil(pd->LP_lower_bound_dual) < (int)ceil(pd->LP_lower_bound))
            ? (int)ceil(pd->LP_lower_bound_dual)
            : (int)ceil(pd->LP_lower_bound);
    pd->LP_lower_bound_BB = CC_MIN(pd->LP_lower_bound, pd->LP_lower_bound_dual);
    pd->LP_lower_min = CC_MIN(pd->LP_lower_min, pd->LP_lower_bound_BB);

    if (pd->iterations % (pd->nb_jobs) == 0 && dbg_lvl() > 0) {
        printf(
            "Current primal LP objective: %19.16f  (LP_dual-bound %19.16f, "
            "lowerbound = %d, eta_in = %f, eta_out = %f).\n",
            pd->LP_lower_bound + pd->off, pd->LP_lower_bound_dual + pd->off,
            pd->lower_bound + pd->off,
            call_get_eta_in(pd->solver_stab) + pd->off,
            pd->LP_lower_bound + pd->off);
    }

CLEAN:
    return val;
}

int solve_relaxation(NodeData* pd) {
    int         val = 0;
    int         status = 0;
    double      real_time_solve_lp = 0.0;
    Statistics* statistics = pd->stat;

    /** Compute LP relaxation */
    real_time_solve_lp = getRealTime();
    CCutil_start_resume_time(&(statistics->tot_solve_lp));
    val = lp_interface_optimize(pd->RMP, &status);
    CCcheck_val_2(val, "lp_interface_optimize failed");
    CCutil_suspend_timer(&(statistics->tot_solve_lp));
    real_time_solve_lp = getRealTime() - real_time_solve_lp;
    statistics->real_time_solve_lp += real_time_solve_lp;

    if (dbg_lvl() > 1) {
        printf("Simplex took %f seconds.\n", real_time_solve_lp);
        fflush(stdout);
    }

    if (dbg_lvl() > 1) {
        print_ages(pd);
    }

    switch (status) {
        case LP_INTERFACE_OPTIMAL:
            /** grow ages of the different columns */
            val = grow_ages(pd);
            CCcheck_val_2(val, "Failed in grow_ages");
            /** get the dual variables and make them feasible */
            val = lp_interface_pi(pd->RMP, &g_array_index(pd->pi, double, 0));
            CCcheck_val_2(val, "lp_interface_pi failed");
            /** Compute the objective function */
            val = compute_objective(pd);
            CCcheck_val_2(val, "Failed in compute_objective");
            break;

        case LP_INTERFACE_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val =
                lp_interface_pi_inf(pd->RMP, &g_array_index(pd->pi, double, 0));
            CCcheck_val_2(val, "Failed at lp_interface_pi_inf");
            break;
    }

CLEAN:

    return val;
}

int compute_lower_bound(NodeData* pd) {
    int         j = 0;
    int         val = 0;
    int         has_cols = 1;
    int         has_cuts = 0;
    int         nb_non_improvements = 0;
    int         status = GRB_LOADED;
    double      real_time_pricing = 0.0;
    Parms*      parms = pd->parms;
    Statistics* statistics = pd->stat;

    if (dbg_lvl() > 1) {
        printf(
            "Starting compute_lower_bound with lb %d and ub %d at depth %d(id "
            "= "
            "%d)\n",
            pd->lower_bound, pd->upper_bound, pd->depth, pd->id);
    }

    CCutil_start_resume_time(&(statistics->tot_lb));

    /**
     * Construction of new solution if localPoolColPool is empty
     */
    if (pd->localColPool->len == 0) {
        // add_solution_to_colpool(problem->opt_sol, pd);
    }
    reset_nb_layers(pd->jobarray);

    if (!pd->RMP) {
        val = build_rmp(pd);
        CCcheck_val(val, "build_lp failed");
    }

    pd->retirementage = (int)sqrt(pd->nb_jobs) + CLEANUP_ITERATION;
    check_schedules(pd);
    delete_infeasible_schedules(pd);

    // solve_relaxation(problem, pd);
    // do {
    do {
        has_cols = 1;
        has_cuts = 0;
        pd->update = 0;
        CCutil_suspend_timer(&(statistics->tot_cputime));
        CCutil_resume_timer(&(statistics->tot_cputime));
        while ((pd->iterations < pd->maxiterations) && has_cols &&
               statistics->tot_cputime.cum_zeit <= parms->branching_cpu_limit) {
            /**
             * Delete old columns
             */
            if (pd->zero_count > pd->nb_jobs * min_nb_del_row_ratio &&
                status == GRB_OPTIMAL) {
                // val = delete_old_schedules(pd);
                CCcheck_val_2(val, "Failed in delete_old_cclasses");
            }
            solve_relaxation(pd);

            /**
             * Solve the pricing problem
             */
            real_time_pricing = getRealTime();
            CCutil_start_resume_time(&statistics->tot_pricing);
            val = lp_interface_status(pd->RMP, &status);
            CCcheck_val_2(val, "Failed in status");

            switch (status) {
                case GRB_OPTIMAL:
                    pd->iterations++;
                    pd->status = infeasible;

                    val = solve_pricing(pd);
                    CCcheck_val_2(val, "Failed in solving pricing");
                    break;

                case GRB_INFEASIBLE:
                    val = solve_farkas_dbl(pd);
                    CCcheck_val_2(val, "Failed in solving farkas");
                    break;
            }

            CCutil_suspend_timer(&statistics->tot_pricing);
            real_time_pricing = getRealTime() - real_time_pricing;
            statistics->real_time_pricing += real_time_pricing;

            switch (status) {
                case GRB_OPTIMAL:
                    has_cols = (call_stopping_criteria(pd->solver_stab) &&
                                (call_get_eta_sep(pd->solver_stab) <
                                 pd->upper_bound - 1.0 + EPS));
                    // pd->nb_new_sets = 0;
                    // || nb_non_improvements > 5;  // ||
                    // (ceil(pd->eta_in - 0.00001) >= pd->eta_out);

                    break;

                case GRB_INFEASIBLE:
                    has_cols = (pd->nb_new_sets == 0);
                    // pd->nb_new_sets = 0;
                    break;
            }

            CCutil_suspend_timer(&(statistics->tot_cputime));
            CCutil_resume_timer(&(statistics->tot_cputime));
        }

        switch (status) {
            case GRB_OPTIMAL:
                /**
                 * change status of problem
                 */
                // if (problem->status == no_sol) {
                //     problem->status = lp_feasible;
                // }

                if (dbg_lvl() > 1) {
                    printf(
                        "Found lb = %d (%f) upper_bound = %d (id= %d, "
                        "iterations = "
                        "%d).\n",
                        pd->lower_bound, pd->LP_lower_bound, pd->upper_bound,
                        pd->id, pd->iterations);
                }

                /**
                 * Compute the objective function
                 */
                // pd->retirementage = 0;
                // delete_old_schedules(pd);
                check_schedules(pd);
                delete_infeasible_schedules(pd);
                solve_relaxation(pd);
                // compute_objective(pd);
                if (pd->localColPool->len > 0) {
                    val = construct_lp_sol_from_rmp(pd);
                    CCcheck_val_2(val, "Failed in construct lp sol from rmp\n");
                    // solve_relaxation(pd);
                    // delete_old_schedules(pd);
                    // delete_unused_rows(pd);
                    // solve_relaxation(pd);
                    // construct_lp_sol_from_rmp(pd);
                    if (!call_is_integer_solution(pd->solver)) {
                        // has_cuts = (generate_cuts(pd) > 0);
                        has_cuts = 0;
                        // call_update_duals(pd->solver_stab);
                        // lp_interface_write(pd->RMP, "test.lp");
                    }
                }
                break;

            case GRB_INFEASIBLE:
                pd->status = infeasible;
                lp_interface_write(pd->RMP, "infeasible_RMP.lp");
                lp_interface_compute_IIS(pd->RMP);
        }
    } while (has_cuts || pd->update);

    if (dbg_lvl() > 0 || pd->id == 0) {
        printf("iterations = %d\n", pd->iterations);
        printf("lowerbound %d\n", pd->lower_bound + pd->off);
        printf("LP value = %f \n", pd->LP_lower_bound + pd->off);
    }

    if (pd->iterations < pd->maxiterations &&
        statistics->tot_cputime.cum_zeit <= parms->branching_cpu_limit) {
    } else {
        switch (status) {
            case GRB_OPTIMAL:
                pd->status = LP_bound_estimated;
                break;

            case GRB_INFEASIBLE:
                pd->status = infeasible;
                break;
        }
    }
    // } while (pd->depth == 1);

    if (pd->depth == 0) {
        statistics->global_lower_bound =
            CC_MAX(pd->lower_bound + pd->off, statistics->global_lower_bound);
        statistics->root_lower_bound = statistics->global_lower_bound;
        statistics->root_upper_bound = statistics->global_upper_bound;
        statistics->root_rel_error =
            (double)(statistics->global_upper_bound -
                     statistics->global_lower_bound) /
            ((double)statistics->global_lower_bound + EPS);
    }

    fflush(stdout);
    statistics->nb_generated_col += pd->iterations;
    CCutil_suspend_timer(&(statistics->tot_lb));

CLEAN:
    return val;
}

int print_x(NodeData* pd) {
    int val = 0;
    int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in lp_interface_status");

    switch (status) {
        case GRB_OPTIMAL:
            val = lp_interface_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            assert(pd->localColPool->len == nb_cols - pd->id_pseudo_schedules);
            pd->lambda = CC_SAFE_REALLOC(
                pd->lambda, nb_cols - pd->id_pseudo_schedules, double);
            CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
            val = lp_interface_x(pd->RMP, pd->lambda, pd->id_pseudo_schedules);
            CCcheck_val_2(val, "Failed in lp_interface_x");

            for (guint i = 0; i < pd->localColPool->len; ++i) {
                GPtrArray* tmp =
                    ((ScheduleSet*)g_ptr_array_index(pd->localColPool, i))
                        ->job_list;
                if (pd->lambda[i] > EPS) {
                    g_ptr_array_foreach(tmp, g_print_machine, NULL);
                    printf("\n");
                }
            }
            break;
    }

CLEAN:

    return val;
}

int calculate_nb_layers(NodeData* pd, int k) {
    int val = 0;
    int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in lp_interface_status");
    if (status == GRB_LOADED) {
        lp_interface_optimize(pd->RMP, &status);
    }

    reset_nb_layers(pd->jobarray);

    switch (status) {
        case GRB_OPTIMAL:
            val = lp_interface_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            assert(pd->localColPool->len == nb_cols - pd->id_pseudo_schedules);
            pd->lambda = CC_SAFE_REALLOC(pd->lambda, nb_cols, double);
            CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
            val = lp_interface_x(pd->RMP, pd->lambda, 0);
            CCcheck_val_2(val, "Failed in lp_interface_x");

            for (unsigned i = 0; i < pd->localColPool->len; ++i) {
                if (pd->lambda[i + pd->id_pseudo_schedules] > EPS) {
                    ScheduleSet* tmp =
                        (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
                    for (int j = 0; j < (int)tmp->job_list->len - k; ++j) {
                        Job* j1 = (Job*)g_ptr_array_index(tmp->job_list, j);
                        Job* j2 = (Job*)g_ptr_array_index(tmp->job_list, j + k);
                        if (j1 == j2) {
                            j1->num_layers = 1;
                            tmp->del = 1;
                        }
                    }
                }
            }
            break;
        default:
            printf("tres test \n");
    }

CLEAN:

    return val;
}

int calculate_x_e(NodeData* pd) {
    int val = 0;
    int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in lp_interface_status");

    switch (status) {
        case GRB_OPTIMAL:
            val = lp_interface_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->lambda = CC_SAFE_REALLOC(pd->lambda, nb_cols, double);
            CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
            val = lp_interface_x(pd->RMP, pd->lambda, 0);
            CCcheck_val_2(val, "Failed in lp_interface_x");
            // pd->x_e =
            //     CC_SAFE_REALLOC(pd->x_e, get_nb_edges(pd->solver), double);
            // CCcheck_NULL_2(pd->x_e, "Failed to reallocate memory to
            // pd->x_e");

            // for (unsigned i = 0; i < get_nb_edges(pd->solver); ++i) {
            //     pd->x_e[i] = 0.0;
            // }
            construct_lp_sol_from_rmp(pd);
            break;
    }

CLEAN:

    return val;
}

int check_schedules(NodeData* pd) {
    int val = 0;
    int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in lp_interface_status");

    val = lp_interface_get_nb_cols(pd->RMP, &nb_cols);
    CCcheck_val_2(val, "Failed to get nb cols");
    assert(nb_cols - pd->id_pseudo_schedules == pd->localColPool->len);
    if (dbg_lvl() > 1) {
        printf("number of cols check %d\n", nb_cols - pd->id_pseudo_schedules);
    }
    for (unsigned i = 0; i < pd->localColPool->len; ++i) {
        ScheduleSet* tmp = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        if (check_schedule_set(tmp, pd) == 1) {
            tmp->del = 1;
        } else {
            tmp->del = 0;
        }
    }

CLEAN:

    return val;
}