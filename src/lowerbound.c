#include <solver.h>
#include <wct.h>
#include "lp.h"
#include "scheduleset.h"
#include "util.h"
#include "wctparms.h"
#include "wctprivate.h"

static const double min_nb_del_row_ratio = 0.9;

void g_print_ages_col(gpointer data, gpointer user_data) {
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

    if (pd->column_status[x->id] == wctlp_LOWER ||
        pd->column_status[x->id] == wctlp_FREE) {
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
    int nb_cols;
    wctlp_get_nb_cols(pd->RMP, &nb_cols);
    assert(nb_cols == pd->localColPool->len);
    CC_IFFREE(pd->column_status, int);
    pd->column_status = (int*)CC_SAFE_MALLOC(nb_cols, int);
    CCcheck_NULL_2(pd->column_status, "Failed to allocate column_status");
    val = wctlp_basis_cols(pd->RMP, pd->column_status, 0);
    CCcheck_val_2(val, "Failed in wctlp_basis_cols");
    pd->zero_count = 0;

    g_ptr_array_foreach(pd->localColPool, g_grow_ages, pd);

CLEAN:
    return val;
}

int delete_old_cclasses(NodeData* pd) {
    int          val = 0;
    int          min_numdel = pd->nb_jobs * min_nb_del_row_ratio;
    int          nb_col;
    guint        i;
    guint        count = pd->localColPool->len;
    ScheduleSet* tmp_schedule;
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
        wctlp_get_nb_cols(pd->RMP, &nb_col);
        assert(nb_col == count);
        for (i = 0; i < count; ++i) {
            tmp_schedule =
                (ScheduleSet*)g_ptr_array_index(pd->localColPool, it);
            if (tmp_schedule->age <= pd->retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = wctlp_deletecols(pd->RMP, first_del, last_del);
                    CCcheck_val_2(val, "Failed in wctlp_deletecols");
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
            wctlp_deletecols(pd->RMP, first_del, last_del);
            CCcheck_val_2(val, "Failed in wctlp_deletecols");
            g_ptr_array_remove_range(pd->localColPool, first_del,
                                     last_del - first_del + 1);
        }

        if (dbg_lvl() > 1) {
            printf("Deleted %d out of %d columns with age > %d.\n",
                   pd->zero_count, count, pd->retirementage);
        }

        wctlp_get_nb_cols(pd->RMP, &nb_col);
        assert(pd->localColPool->len == nb_col);
        for (i = 0; i < pd->localColPool->len; ++i) {
            tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
            tmp_schedule->id = i;
        }
        pd->zero_count = 0;
    }

CLEAN:
    return val;
}

int delete_infeasible_cclasses(NodeData* pd) {
    int          val = 0;
    int          nb_col;
    guint        i;
    guint        count = pd->localColPool->len;
    ScheduleSet* tmp_schedule;
    /** pd->zero_count can be deprecated! */
    pd->zero_count = 0;

    for (i = pd->nb_jobs; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        if (tmp_schedule->age > 0) {
            pd->zero_count++;
        }
    }


    int it = 0;
    int first_del = pd->nb_jobs-1;
    int last_del = pd->nb_jobs-1;
    for (i = pd->nb_jobs; i < count; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, it);
        if (tmp_schedule->del != 1) {
            if (first_del != pd->nb_jobs-1) {
                /** Delete recently found deletion range.*/
                val = wctlp_deletecols(pd->RMP, first_del, last_del);
                CCcheck_val_2(val, "Failed in wctlp_deletecols");
                g_ptr_array_remove_range(pd->localColPool, first_del,
                                         last_del - first_del + 1);
                it = it - (last_del - first_del);
                first_del = last_del = pd->nb_jobs-1;
            } else {
                it++;
            }
        } else {
            if (first_del == pd->nb_jobs-1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            it++;
        }
    }

    if (first_del != pd->nb_jobs-1) {
        wctlp_deletecols(pd->RMP, first_del, last_del);
        CCcheck_val_2(val, "Failed in wctlp_deletecols");
        g_ptr_array_remove_range(pd->localColPool, first_del,
                                 last_del - first_del + 1);
    }

    if (dbg_lvl() > 1) {
        printf("Deleted %d out of %d columns with age > %d.\n", pd->zero_count,
               count, pd->retirementage);
    }

    wctlp_get_nb_cols(pd->RMP, &nb_col);
    assert(pd->localColPool->len == nb_col);
    printf("number of cols = %d\n", nb_col);
    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        tmp_schedule->id = i;
    }
    printf("tset");
    pd->zero_count = 0;

CLEAN:
    return val;
}

void g_make_pi_feasible(gpointer data, gpointer user_data) {
    ScheduleSet* x = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;
    Job*         tmp_j;

    int    i;
    double colsum = .0;

    for (i = 0; i < x->job_list->len; ++i) {
        tmp_j = (Job*)g_ptr_array_index(x->job_list, i);
        if (signbit((double)pd->pi->data[tmp_j->job])) {
            pd->pi->data[tmp_j->job]= 0.0;
        }

        colsum += pd->pi->data[tmp_j->job];
        colsum = nextafter(colsum, DBL_MAX);
    }

    if (!signbit((double)pd->pi->data[pd->nb_jobs])) {
        g_array_index(pd->pi, double, pd->nb_jobs)  = 0.0;
    }

    colsum += (double)pd->pi->data[pd->nb_jobs];
    colsum = nextafter(colsum, DBL_MAX);

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (i = 0; i < x->job_list->len; ++i) {
            tmp_j = (Job*)g_ptr_array_index(x->job_list, i);
            g_array_index(pd->pi, double, tmp_j->job)  /= colsum;
            g_array_index(pd->pi, double, tmp_j->job) *= x->total_weighted_completion_time;
            newcolsum += g_array_index(pd->pi, double, tmp_j->job) ;
        }

        g_array_index(pd->pi, double, pd->nb_jobs) /= colsum;
        g_array_index(pd->pi, double, pd->nb_jobs)  *= x->total_weighted_completion_time;
        newcolsum += g_array_index(pd->pi, double, pd->nb_jobs) ;

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
    Job*         tmp_j;

    int    i;
    double colsum = .0;

    for (i = 0; i < x->job_list->len; ++i) {
        double *tmp = &g_array_index(pd->pi, double, i);
        tmp_j = (Job*)g_ptr_array_index(x->job_list, i);
        if (signbit(*tmp)) {
            *tmp = 0.0;
        }

        colsum += *tmp;
        colsum = nextafter(colsum, DBL_MAX);
    }

    colsum +=  g_array_index(pd->pi,double, pd->nb_jobs) ;

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (i = 0; i < x->job_list->len; ++i) {
            tmp_j = (Job*)g_ptr_array_index(x->job_list, i);
            double *tmp = &g_array_index(pd->pi, double, i);
            *tmp /= colsum;
            *tmp *= x->total_weighted_completion_time;
            newcolsum += *tmp;
        }

        double *tmp = &g_array_index(pd->pi, double, pd->nb_jobs);
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

int compute_objective(NodeData* pd, Parms* parms) {
    int val = 0;
    int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    double *tmp = &g_array_index(pd->pi, double, 0);
    double *tmp_rhs = &g_array_index(pd->rhs, double, 0);
    for (i = 0; i < pd->nb_jobs + 1; i++) {
        pd->LP_lower_bound_dual += tmp[i] * tmp_rhs[i];
    }
    pd->LP_lower_bound_dual -= 0.0001;

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(pd->RMP, &(pd->LP_lower_bound));
    pd->LP_lower_bound -= 0.0001;
    CCcheck_val_2(val, "wctlp_objval failed");
    pd->lower_bound =
        ((int)ceil(pd->LP_lower_bound_dual) < (int)ceil(pd->LP_lower_bound))
            ? (int)ceil(pd->LP_lower_bound_dual)
            : (int)ceil(pd->LP_lower_bound);
    pd->LP_lower_bound_BB = CC_MIN(pd->LP_lower_bound, pd->LP_lower_bound_dual);
    pd->eta_out = pd->LP_lower_bound_BB;

    if (pd->iterations % pd->nb_jobs == 0) {
        printf(
            "Current primal LP objective: %19.16f  (LP_dual-bound %19.16f, "
            "lowerbound = %d, eta_in = %f, eta_out = %f).\n",
            pd->LP_lower_bound, pd->LP_lower_bound_dual, pd->lower_bound,
            pd->eta_in, pd->eta_out);
    }

CLEAN:
    return val;
}

int solve_relaxation(Problem* problem, NodeData *pd) {
    int val = 0;
    int status;
    Parms *parms = &(problem->parms);
    double real_time_solve_lp;

    /** Compjute LP relaxation */
    real_time_solve_lp = getRealTime();
    CCutil_start_resume_time(&(problem->tot_solve_lp));
    val = wctlp_optimize(pd->RMP, &status);
    CCcheck_val_2(val, "wctlp_optimize failed");
    CCutil_suspend_timer(&(problem->tot_solve_lp));
    real_time_solve_lp = getRealTime() - real_time_solve_lp;
    problem->real_time_solve_lp += real_time_solve_lp;

    if (dbg_lvl() > 1) {
        printf("Simplex took %f seconds.\n", real_time_solve_lp);
        fflush(stdout);
    }

    if (dbg_lvl() > 1) {
        print_ages(pd);
    }

    switch (status) {
        case WCTLP_OPTIMAL:
            /** grow ages of the different columns */
            val = grow_ages(pd);
            CCcheck_val_2(val, "Failed in grow_ages");
            /** get the dual variables and make them feasible */
            val = wctlp_pi(pd->RMP, &g_array_index(pd->pi,double, 0));
            CCcheck_val_2(val, "wctlp_pi failed");
            /** Compute the objective function */
            val = compute_objective(pd, parms);
            CCcheck_val_2(val, "Failed in compute_objective");
            memcpy(&g_array_index(pd->pi_out,double,0), &g_array_index(pd->pi, double, 0), sizeof(double) * (pd->nb_jobs + 1));
            pd->eta_out = pd->LP_lower_bound_dual;
            break;

        case WCTLP_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = wctlp_pi_inf(pd->RMP, &g_array_index(pd->pi, double, 0));
            CCcheck_val_2(val, "Failed at wctlp_pi_inf");
            break;
    }

    CLEAN:

    return val;
}

int compute_lower_bound(Problem* problem, NodeData* pd) {
    int    j, val = 0;
    int    break_while_loop = 1;
    int    nb_non_improvements = 0;
    int    status = GRB_LOADED;
    double real_time_pricing;
    Parms* parms = &(problem->parms);

    if (dbg_lvl() > 1) {
        printf(
            "Starting compute_lower_bound with lb %d and ub %d at depth %d(id "
            "= "
            "%d, opt_track = %d)\n",
            pd->lower_bound, pd->upper_bound, pd->depth, pd->id, pd->opt_track);
    }

    CCutil_start_resume_time(&(problem->tot_lb));

    /**
     * Construction of new solution if localPoolColPool is empty
     */
    if (pd->localColPool->len == 0) {
        // add_solution_to_colpool(problem->opt_sol, pd);
    }
    reset_nb_layers(pd->jobarray);

    if (!pd->RMP) {
        val = build_rmp(pd, 0);
        CCcheck_val(val, "build_lp failed");
    }

    pd->retirementage = (int)sqrt(pd->nb_jobs) + 30;
    // check_schedules(pd);
    // delete_infeasible_cclasses(pd);

    /** Init alpha */
    switch (parms->stab_technique) {
        case stab_wentgnes:
            pd->alpha = parms->alpha;
            break;

        case stab_dynamic:
            pd->alpha = 0.0;
            break;

        case no_stab:
            break;
    }

    // solve_relaxation(problem, pd);

    break_while_loop = 0;
    CCutil_suspend_timer(&(problem->tot_cputime));
    CCutil_resume_timer(&(problem->tot_cputime));

    while ((pd->iterations < pd->maxiterations) && !break_while_loop &&
           problem->tot_cputime.cum_zeit <=
               problem->parms.branching_cpu_limit) {
        /**
         * Delete old columns
         */
        // if (pd->zero_count > pd->nb_jobs * min_nb_del_row_ratio &&
        //     status == GRB_OPTIMAL) {
        //     val = delete_old_cclasses(pd);
        //     CCcheck_val_2(val, "Failed in delete_old_cclasses");
        // }

        /**
         * Solve the pricing problem
         */
        real_time_pricing = getRealTime();
        CCutil_start_resume_time(&problem->tot_pricing);
        val = wctlp_status(pd->RMP, &status);
        CCcheck_val_2(val, "Failed in status");

        switch (status) {
            case GRB_OPTIMAL:
                pd->iterations++;
                pd->status = infeasible;

                if (pd->iterations < pd->maxiterations) {
                    switch (parms->stab_technique) {
                        case stab_wentgnes:
                            val = solve_stab(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case stab_dynamic:
                            val = solve_stab_dynamic(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case stab_hybrid:
                            val = solve_stab_hybrid(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab_hybrid");
                            break;

                        case no_stab:
                            val = solve_pricing(pd, parms, 0);
                            CCcheck_val_2(val, "Failed in solving pricing");
                            break;
                    }
                }

                break;

            case GRB_INFEASIBLE:
                val = solve_farkas_dbl(pd);
                CCcheck_val_2(val, "Failed in solving farkas");
                break;
        }

        CCutil_suspend_timer(&problem->tot_pricing);
        real_time_pricing = getRealTime() - real_time_pricing;
        problem->real_time_pricing += real_time_pricing;

        if (pd->update) {
            for (j = 0; j < pd->nb_new_sets; j++) {
                val = add_lhs_scheduleset_to_rmp(pd->newsets + j, pd);
                CCcheck_val_2(val, "wctlp_addcol failed");
                g_ptr_array_add(pd->localColPool, pd->newsets + j);
            }
            pd->newsets = NULL;
            nb_non_improvements = 0;
        } else {
            nb_non_improvements++;
        }

        switch (status) {
            case GRB_OPTIMAL:
                switch (parms->stab_technique) {
                    case stab_wentgnes:
                    case stab_dynamic:
                    case stab_hybrid:
                        break_while_loop =
                            (CC_ABS(pd->eta_out - pd->eta_in) < 0.00001) ||
                            nb_non_improvements > 5;  // || (ceil(pd->eta_in - 0.00001) >=
                                    // pd->eta_out);
                        break;

                    case no_stab:
                        break_while_loop =
                            (pd->nb_new_sets == 0 || nb_non_improvements > 5);
                        break;
                }

                break;

            case GRB_INFEASIBLE:
                break_while_loop = (pd->nb_new_sets == 0);
                break;
        }

        solve_relaxation(problem, pd);

        CCutil_suspend_timer(&(problem->tot_cputime));
        CCutil_resume_timer(&(problem->tot_cputime));
    }

    if (pd->iterations < pd->maxiterations &&
        problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        switch (status) {
            case GRB_OPTIMAL:
                /**
                 * change status of problem
                 */
                if (problem->status == no_sol) {
                    problem->status = lp_feasible;
                }

                if (dbg_lvl() > 1) {
                    printf(
                        "Found lb = %d (%f) upper_bound = %d (id= %d, "
                        "iterations = "
                        "%d,opt_track = %d).\n",
                        pd->lower_bound, pd->LP_lower_bound, pd->upper_bound,
                        pd->id, pd->iterations, pd->opt_track);
                }

                if(parms->reduce_cost_fixing == yes_reduced_cost) {
                    CCutil_start_resume_time(&(problem->tot_reduce_cost_fixing));
                    reduce_cost_fixing(pd);
                    CCutil_suspend_timer(&(problem->tot_reduce_cost_fixing));
                    // print_interval_pair(pd->ordered_jobs);
                }

                /**
                 * Compute the objective function
                 */
                val = wctlp_optimize(pd->RMP, &status);
                CCcheck_val_2(val, "wctlp_optimize failed");
                val = compute_objective(pd, parms);
                CCcheck_val_2(val, "Failed in compute_objective");
                // construct_lp_sol_from_rmp(pd);
                memcpy(&g_array_index(pd->pi_out,double,0), &g_array_index(pd->pi, double, 0) , sizeof(double) * (pd->nb_jobs + 1));
                printf("size evolution %lu\n", get_nb_vertices(pd->solver));
                break;

            case GRB_INFEASIBLE:
                pd->status = infeasible;
                pd->test = 0;
                wctlp_write(pd->RMP, "infeasible_RMP.lp");
                wctlp_compute_IIS(pd->RMP);
        }
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

    if (dbg_lvl() > -1) {
        printf("iterations = %d\n", pd->iterations);
        printf("lowerbound %d\n", pd->lower_bound + pd->problem->off);
        printf("LP value = %f \n", pd->LP_lower_bound + pd->problem->off);
    }
    problem->global_lower_bound =
        CC_MAX(pd->lower_bound + pd->problem->off, problem->global_lower_bound);

    if (pd == &(problem->root_pd)) {
        problem->root_lower_bound = problem->global_lower_bound;
        problem->root_upper_bound = problem->global_upper_bound;
        problem->root_rel_error =
            (double)(problem->global_upper_bound -
                     problem->global_lower_bound) /
            ((double)problem->global_lower_bound + 0.000001);
    }

    fflush(stdout);
    problem->nb_generated_col += pd->iterations;
    CCutil_suspend_timer(&(problem->tot_lb));

CLEAN:
    return val;
}

int print_x(NodeData* pd) {
    int val = 0;
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status");

    switch (status) {
        case GRB_OPTIMAL:
            val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->lambda = CC_SAFE_REALLOC(pd->lambda, nb_cols, double);
            CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
            val = wctlp_x(pd->RMP, pd->lambda, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");

            for (int i = 0; i < nb_cols; ++i) {
                GPtrArray* tmp =
                    ((ScheduleSet*)g_ptr_array_index(pd->localColPool, i))
                        ->job_list;
                if (pd->lambda[i] > 0.00001) {
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
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status");
    if (status == GRB_LOADED) {
        wctlp_optimize(pd->RMP, &status);
    }

    reset_nb_layers(pd->jobarray);

    switch (status) {
        case GRB_OPTIMAL:
            val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->lambda = CC_SAFE_REALLOC(pd->lambda, nb_cols, double);
            CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
            val = wctlp_x(pd->RMP, pd->lambda, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");

            for (unsigned i = 0; i < nb_cols; ++i) {
                if (pd->lambda[i] > 0.00001) {
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
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status")

        switch (status) {
        case GRB_OPTIMAL:
            val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->lambda = CC_SAFE_REALLOC(pd->lambda, nb_cols, double);
            CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
            val = wctlp_x(pd->RMP, pd->lambda, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");
            pd->x_e =
                CC_SAFE_REALLOC(pd->x_e, get_nb_edges(pd->solver), double);
            CCcheck_NULL_2(pd->x_e, "Failed to reallocate memory to  pd->x_e");

            for (unsigned i = 0; i < get_nb_edges(pd->solver); ++i) {
                pd->x_e[i] = 0.0;
            }
            construct_lp_sol_from_rmp(pd);
            break;
    }

CLEAN:

    return val;
}

int check_schedules(NodeData* pd) {
    int val = 0;
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status")

    val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
    CCcheck_val_2(val, "Failed to get nb cols");
    assert(nb_cols == pd->localColPool->len);
    printf("number of cols check %d\n", nb_cols);
    for (unsigned i = 0; i < nb_cols; ++i) {
        ScheduleSet* tmp = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        if (check_schedule_set(tmp, pd) == 1) {
            tmp->del = 0;
        } else {
            tmp->del = 1;
            // printf("deleted column\n");
            // print_schedule(tmp, 1);
        }
    }


CLEAN:

    return val;
}
