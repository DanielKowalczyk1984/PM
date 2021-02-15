#include "gurobi_c.h"
#include "job.h"
#include "lp.h"
#include "parms.h"
#include "scheduleset.h"
#include "solver.h"
#include "util.h"
#include "wctprivate.h"

static const double min_nb_del_row_ratio = 0.9;

void g_print_ages_col(gpointer data, MAYBE_UNUSED gpointer user_data) {
    ScheduleSet* x = static_cast<ScheduleSet*>(data);

    fmt::print(" {}", x->age);
}

/** Help function for column generation */
void NodeData::print_ages() {
    fmt::print("AGES:");

    g_ptr_array_foreach(localColPool, g_print_ages_col, NULL);

    fmt::print("\n");
}

void g_grow_ages(gpointer data, gpointer user_data) {
    ScheduleSet* x = static_cast<ScheduleSet*>(data);
    NodeData*    pd = static_cast<NodeData*>(user_data);

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

int NodeData::grow_ages() {
    int val = 0;
    int nb_cols = 0;
    lp_interface_get_nb_cols(RMP, &nb_cols);
    assert(nb_cols - id_pseudo_schedules == localColPool->len);
    CC_IFFREE(column_status, int);
    if (localColPool->len > 0) {
        column_status = CC_SAFE_MALLOC(localColPool->len, int);
        CCcheck_NULL_2(column_status, "Failed to allocate column_status");
        val = lp_interface_basis_cols(RMP, column_status, id_pseudo_schedules);
        CCcheck_val_2(val, "Failed in lp_interface_basis_cols");
        zero_count = 0;

        g_ptr_array_foreach(localColPool, g_grow_ages, this);
    }

CLEAN:
    return val;
}

int NodeData::delete_unused_rows() {
    int     val = 0;
    int     nb_rows = 0;
    double* slack_tmp = &g_array_index(slack, double, 0);
    GArray* del_indices = g_array_new(FALSE, FALSE, sizeof(int));

    lp_interface_get_nb_rows(RMP, &nb_rows);
    // assert(nb_rows == pd->slack->len);
    lp_interface_slack(RMP, slack_tmp);

    int it = id_valid_cuts;
    int first_del = -1;
    int last_del = -1;
    for (int i = id_valid_cuts; i < nb_rows; i++) {
        if (std::fabs(slack_tmp[i]) < EPS) {
            if (first_del != -1) {
                val = delete_unused_rows_range(first_del, last_del);
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
        delete_unused_rows_range(first_del, last_del);
    }

    call_update_rows_coeff();

    return val;
}

void g_scheduleset_count_zero(gpointer data, gpointer user_data) {
    ScheduleSet* tmp = static_cast<ScheduleSet*>(data);
    int*         aux = static_cast<int*>(user_data);

    if (tmp->age > 0) {
        (*aux)++;
    }
}

int NodeData::delete_old_schedules() {
    int          val = 0;
    int          min_numdel = floor(nb_jobs * min_nb_del_row_ratio);
    int          nb_col = 0;
    guint        i = 0U;
    guint        count = localColPool->len;
    ScheduleSet* tmp_schedule = (ScheduleSet*)NULL;
    /** pd->zero_count can be deprecated! */
    zero_count = 0;

    g_ptr_array_foreach(localColPool, g_scheduleset_count_zero, &zero_count);

    if (zero_count > min_numdel) {
        int it = 0;
        int first_del = -1;
        int last_del = -1;
        lp_interface_get_nb_cols(RMP, &nb_col);
        assert(nb_col - id_pseudo_schedules == count);
        for (i = 0; i < count; ++i) {
            tmp_schedule =
                static_cast<ScheduleSet*>(g_ptr_array_index(localColPool, it));
            if (tmp_schedule->age <= retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = lp_interface_deletecols(
                        RMP, first_del + id_pseudo_schedules,
                        last_del + id_pseudo_schedules);
                    // CCcheck_val_2(val, "Failed in lp_interface_deletecols");
                    g_ptr_array_remove_range(localColPool, first_del,
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
            lp_interface_deletecols(RMP, first_del + id_pseudo_schedules,
                                    last_del + id_pseudo_schedules);
            g_ptr_array_remove_range(localColPool, first_del,
                                     last_del - first_del + 1);
        }

        if (dbg_lvl() > 1) {
            fmt::print("Deleted %d out of %d columns with age > %d.\n",
                       zero_count, count, retirementage);
        }

        lp_interface_get_nb_cols(RMP, &nb_col);
        assert(localColPool->len == nb_col - id_pseudo_schedules);
        for (i = 0; i < localColPool->len; ++i) {
            tmp_schedule = (ScheduleSet*)g_ptr_array_index(localColPool, i);
            tmp_schedule->id = i;
        }
        zero_count = 0;
    }

CLEAN:
    return val;
}

int NodeData::delete_infeasible_schedules() {
    int          val = 0;
    int          nb_col = 0;
    guint        i = 0;
    guint        count = localColPool->len;
    ScheduleSet* tmp_schedule = (ScheduleSet*)NULL;
    /** pd->zero_count can be deprecated! */
    zero_count = 0;

    int it = 0;
    int first_del = -1;
    int last_del = -1;
    lp_interface_get_nb_cols(RMP, &nb_col);
    assert(nb_col - id_pseudo_schedules == count);
    for (i = 0; i < count; ++i) {
        // while (it < localColPool->len) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(localColPool, it);
        if (tmp_schedule->del != 0) {
            if (first_del != -1) {
                /** Delete recently found deletion range.*/
                val = lp_interface_deletecols(RMP,
                                              first_del + id_pseudo_schedules,
                                              last_del + id_pseudo_schedules);
                CCcheck_val_2(val, "Failed in lp_interface_deletecols");
                g_ptr_array_remove_range(localColPool, first_del,
                                         last_del - first_del + 1);
                zero_count += last_del - first_del + 1;
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
        lp_interface_deletecols(RMP, first_del + id_pseudo_schedules,
                                last_del + id_pseudo_schedules);
        CCcheck_val_2(val, "Failed in lp_interface_deletecols");
        g_ptr_array_remove_range(localColPool, first_del,
                                 last_del - first_del + 1);
        zero_count += last_del - first_del + 1;
    }

    if (dbg_lvl() > 1) {
        fmt::print(
            "Deleted {} out of {} columns(infeasible columns after reduce cost "
            "fixing).\n",
            zero_count, count);
    }

    lp_interface_get_nb_cols(RMP, &nb_col);
    assert(localColPool->len == nb_col - id_pseudo_schedules);
    if (dbg_lvl() > 1) {
        fmt::print("number of cols = {}\n", nb_col - id_pseudo_schedules);
    }

    for (i = 0; i < localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(localColPool, i);
        tmp_schedule->id = i;
    }

    if (zero_count > 0) {
        solve_relaxation();
        update = 1;
    }

CLEAN:
    return val;
}

void g_make_pi_feasible(gpointer data, gpointer user_data) {
    ScheduleSet* x = static_cast<ScheduleSet*>(data);
    NodeData*    pd = static_cast<NodeData*>(user_data);
    Job*         tmp_j = nullptr;

    double colsum = .0;

    for (guint i = 0; i < x->job_list->len; ++i) {
        tmp_j = static_cast<Job*>(g_ptr_array_index(x->job_list, i));
        if (signbit((double)pd->pi->data[tmp_j->job])) {
            pd->pi->data[tmp_j->job] = 0.0;
        }

        colsum += pd->pi->data[tmp_j->job];
        colsum = nextafter(colsum, DBL_MAX);
    }

    if (!signbit((double)pd->pi->data[pd->nb_jobs])) {
        g_array_index(pd->pi, double, pd->nb_jobs) = 0.0;
    }

    colsum += static_cast<double>(pd->pi->data[pd->nb_jobs]);
    colsum = nextafter(colsum, DBL_MAX);

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (guint i = 0; i < x->job_list->len; ++i) {
            tmp_j = static_cast<Job*>(g_ptr_array_index(x->job_list, i));
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
            fmt::print(
                R"(Decreased column sum of {} from  {:30.20f} to  {:30.20f}
)",
                x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible(NodeData* pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible, pd);
}

void g_make_pi_feasible_farkas(gpointer data, gpointer user_data) {
    ScheduleSet* x = static_cast<ScheduleSet*>(data);
    NodeData*    pd = static_cast<NodeData*>(user_data);

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
            fmt::print(
                R"(Decreased column sum of {} from  {:30.20f} to  {:30.20f}
)",
                x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void NodeData::make_pi_feasible_farkas_pricing() {
    g_ptr_array_foreach(localColPool, g_make_pi_feasible_farkas, this);
}

int NodeData::compute_objective() {
    int val = 0;
    // int i;
    LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    double* tmp = &g_array_index(pi, double, 0);
    double* tmp_rhs = &g_array_index(rhs, double, 0);

    for (int i = 0; i < nb_rows; i++) {
        // if (i != nb_jobs) {
        // eta_out += tmp[i] * tmp_rhs[i];
        // }
        LP_lower_bound_dual += tmp[i] * tmp_rhs[i];
    }
    LP_lower_bound_dual -= EPS_BOUND;

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = lp_interface_objval(RMP, &(LP_lower_bound));
    LP_lower_bound -= EPS_BOUND;
    // CCcheck_val_2(val, "lp_interface_objval failed");
    lower_bound = (ceil(LP_lower_bound_dual) < ceil(LP_lower_bound))
                      ? static_cast<int>(ceil(LP_lower_bound_dual))
                      : static_cast<int>(ceil(LP_lower_bound));
    LP_lower_bound_BB = std::min(LP_lower_bound, LP_lower_bound_dual);
    LP_lower_min = std::min(LP_lower_min, LP_lower_bound_BB);

    if (iterations % (nb_jobs) == 0 && dbg_lvl() > 0) {
        fmt::print(
            "Current primal LP objective: {:19.16f}  (LP_dual-bound {:19.16f}, "
            "lowerbound = {}, eta_in = {}, eta_out = {}).\n",
            LP_lower_bound + off, LP_lower_bound_dual + off, lower_bound + off,
            solver_stab->get_eta_in() + off, LP_lower_bound + off);
    }

CLEAN:
    return val;
}

int NodeData::solve_relaxation() {
    int         val = 0;
    int         status = 0;
    double      real_time_solve_lp = 0.0;
    Statistics* statistics = stat;

    /** Compute LP relaxation */
    real_time_solve_lp = getRealTime();
    CCutil_start_resume_time(&(statistics->tot_solve_lp));
    val = lp_interface_optimize(RMP, &status);
    // CCcheck_val_2(val, "lp_interface_optimize failed");
    CCutil_suspend_timer(&(statistics->tot_solve_lp));
    real_time_solve_lp = getRealTime() - real_time_solve_lp;
    statistics->real_time_solve_lp += real_time_solve_lp;

    if (dbg_lvl() > 1) {
        fmt::print("Simplex took {} seconds.\n", real_time_solve_lp);
        fflush(stdout);
    }

    if (dbg_lvl() > 1) {
        print_ages();
    }

    switch (status) {
        case LP_INTERFACE_OPTIMAL:
            /** grow ages of the different columns */
            val = grow_ages();
            // CCcheck_val_2(val, "Failed in grow_ages");
            /** get the dual variables and make them feasible */
            val = lp_interface_pi(RMP, &g_array_index(pi, double, 0));
            // CCcheck_val_2(val, "lp_interface_pi failed");
            /** Compute the objective function */
            val = compute_objective();
            // CCcheck_val_2(val, "Failed in compute_objective");
            break;

        case LP_INTERFACE_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = lp_interface_pi_inf(RMP, &g_array_index(pi, double, 0));
            // CCcheck_val_2(val, "Failed at lp_interface_pi_inf");
            break;
    }

CLEAN:

    return val;
}

int NodeData::compute_lower_bound() {
    int         j = 0;
    int         val = 0;
    int         has_cols = 1;
    int         has_cuts = 0;
    int         nb_non_improvements = 0;
    int         status_RMP = GRB_LOADED;
    double      real_time_pricing = 0.0;
    Statistics* statistics = stat;

    if (dbg_lvl() > 1) {
        fmt::print(
            R"(Starting compute_lower_bound with lb {} and ub %d at depth {}(id = {})
)",
            lower_bound, upper_bound, depth, id);
    }

    CCutil_start_resume_time(&(statistics->tot_lb));

    /**
     * Construction of new solution if localPoolColPool is empty
     */
    if (localColPool->len == 0) {
        // add_solution_to_colpool(problem->opt_sol, pd);
    }
    reset_nb_layers(jobarray);

    if (!RMP) {
        val = build_rmp();
        // CCcheck_val(val, "build_lp failed");
    }

    retirementage = static_cast<int>(sqrt(nb_jobs)) + CLEANUP_ITERATION;
    check_schedules();
    delete_infeasible_schedules();

    // solve_relaxation(problem, pd);
    // do {
    do {
        has_cols = 1;
        has_cuts = 0;
        update = 0;
        CCutil_suspend_timer(&(statistics->tot_cputime));
        CCutil_resume_timer(&(statistics->tot_cputime));
        while ((iterations < maxiterations) && has_cols &&
               statistics->tot_cputime.cum_zeit <= parms->branching_cpu_limit) {
            /**
             * Delete old columns
             */
            if (zero_count > nb_jobs * min_nb_del_row_ratio &&
                status == GRB_OPTIMAL) {
                // val = delete_old_schedules(pd);
                // CCcheck_val_2(val, "Failed in delete_old_cclasses");
            }
            solve_relaxation();

            /**
             * Solve the pricing problem
             */
            real_time_pricing = getRealTime();
            CCutil_start_resume_time(&statistics->tot_pricing);
            val = lp_interface_status(RMP, &status_RMP);
            // CCcheck_val_2(val, "Failed in status");

            switch (status_RMP) {
                case GRB_OPTIMAL:
                    iterations++;
                    status = infeasible;

                    val = solve_pricing();
                    // CCcheck_val_2(val, "Failed in solving pricing");
                    break;

                case GRB_INFEASIBLE:
                    val = solve_farkas_dbl();
                    // CCcheck_val_2(val, "Failed in solving farkas");
                    break;
            }

            CCutil_suspend_timer(&statistics->tot_pricing);
            real_time_pricing = getRealTime() - real_time_pricing;
            statistics->real_time_pricing += real_time_pricing;

            switch (status_RMP) {
                case GRB_OPTIMAL:
                    has_cols = (solver_stab->stopping_criteria() &&
                                solver_stab->get_eta_in() <
                                    upper_bound - 1.0 + EPS_BOUND);
                    // nb_new_sets = 0;
                    // || nb_non_improvements > 5;  // ||
                    // (ceil(eta_in - 0.00001) >= eta_out);

                    break;

                case GRB_INFEASIBLE:
                    has_cols = (nb_new_sets == 0);
                    // nb_new_sets = 0;
                    break;
            }

            CCutil_suspend_timer(&(statistics->tot_cputime));
            CCutil_resume_timer(&(statistics->tot_cputime));
        }

        switch (status_RMP) {
            case GRB_OPTIMAL:
                /**
                 * change status of problem
                 */
                // if (problem->status == no_sol) {
                //     problem->status = lp_feasible;
                // }

                if (dbg_lvl() > 1) {
                    fmt::print(
                        R"(Found lb = {} ({}) upper_bound = {} (id = {}, iterations = {}).
)",
                        lower_bound, LP_lower_bound, upper_bound, id,
                        iterations);
                }

                /**
                 * Compute the objective function
                 */
                // retirementage = 0;
                // delete_old_schedules(pd);
                solve_relaxation();
                // double obj;
                // lp_interface_objval(pd->RMP, &obj);
                // printf("test objval = %f\n", obj);
                // check_schedules(pd);
                // delete_infeasible_schedules(pd);
                // solve_relaxation(pd);
                // lp_interface_objval(pd->RMP, &obj);
                // printf("test objval = %f\n", obj);
                // printf("----------------\n");
                // compute_objective(pd);
                if (localColPool->len > 0) {
                    val = construct_lp_sol_from_rmp();
                    // CCcheck_val_2(val, "Failed in construct lp sol from
                    // rmp\n"); solve_relaxation(pd); delete_old_schedules(pd);
                    // delete_unused_rows(pd);
                    // solve_relaxation(pd);
                    // construct_lp_sol_from_rmp(pd);
                    // if (!call_is_integer_solution(pd->solver)) {
                    // has_cuts = (generate_cuts(pd) > 0);
                    // has_cuts = 0;
                    // call_update_duals(pd->solver_stab);
                    // lp_interface_write(pd->RMP, "test.lp");
                    // }
                }
                break;

            case GRB_INFEASIBLE:
                status = infeasible;
                lp_interface_write(RMP, "infeasible_RMP.lp");
                lp_interface_compute_IIS(RMP);
        }
    } while (0);

    if (iterations < maxiterations &&
        statistics->tot_cputime.cum_zeit <= parms->branching_cpu_limit) {
    } else {
        switch (status_RMP) {
            case GRB_OPTIMAL:
                status = LP_bound_estimated;
                break;

            case GRB_INFEASIBLE:
                status = infeasible;
                break;
        }
    }
    // } while (depth == 1);

    if (depth == 0) {
        statistics->global_lower_bound =
            CC_MAX(lower_bound + off, statistics->global_lower_bound);
        statistics->root_lower_bound = statistics->global_lower_bound;
        statistics->root_upper_bound = statistics->global_upper_bound;
        statistics->root_rel_error =
            static_cast<double>(statistics->global_upper_bound -
                                statistics->global_lower_bound) /
            (statistics->global_lower_bound + EPS);
    }

    fflush(stdout);
    statistics->nb_generated_col += iterations;
    CCutil_suspend_timer(&(statistics->tot_lb));

CLEAN:
    return val;
}

int NodeData::print_x() {
    int val = 0;
    int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(RMP, &status);
    CCcheck_val_2(val, "Failed in lp_interface_status");

    switch (status) {
        case GRB_OPTIMAL:
            val = lp_interface_get_nb_cols(RMP, &nb_cols);
            // CCcheck_val_2(val, "Failed to get nb cols");
            assert(localColPool->len == nb_cols - id_pseudo_schedules);
            lambda =
                CC_SAFE_REALLOC(lambda, nb_cols - id_pseudo_schedules, double);
            // CCcheck_NULL_2(lambda, "Failed to allocate memory to x");
            val = lp_interface_x(RMP, lambda, id_pseudo_schedules);
            // CCcheck_val_2(val, "Failed in lp_interface_x");

            for (guint i = 0; i < localColPool->len; ++i) {
                GPtrArray* tmp =
                    ((ScheduleSet*)g_ptr_array_index(localColPool, i))
                        ->job_list;
                if (lambda[i] > EPS) {
                    g_ptr_array_foreach(tmp, g_print_machine, NULL);
                    fmt::print("\n");
                }
            }
            break;
    }

CLEAN:

    return val;
}

int NodeData::check_schedules() {
    int val = 0;
    // int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(RMP, &status);
    // CCcheck_val_2(val, "Failed in lp_interface_status");

    val = lp_interface_get_nb_cols(RMP, &nb_cols);
    // CCcheck_val_2(val, "Failed to get nb cols");
    assert(nb_cols - id_pseudo_schedules == localColPool->len);
    if (dbg_lvl() > 1) {
        fmt::print("number of cols check {}\n", nb_cols - id_pseudo_schedules);
    }
    for (unsigned i = 0; i < localColPool->len; ++i) {
        ScheduleSet* tmp = (ScheduleSet*)g_ptr_array_index(localColPool, i);
        if (check_schedule_set(tmp) == 1) {
            tmp->del = 1;
        } else {
            tmp->del = 0;
        }
    }

    return val;
}