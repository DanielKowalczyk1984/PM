#include <wct.h>
#include <solver.h>

static const double min_ndelrow_ratio = 0.9;

void g_print_ages_col(gpointer data, gpointer user_data) {
    ScheduleSet *x = (ScheduleSet *)data;

    printf(" %4d", x->age);
}

/** Help function for column generation */
static void print_ages(NodeData *pd) {
    printf("AGES:");

    g_ptr_array_foreach(pd->localColPool, g_print_ages_col, NULL);

    printf("\n");
}

void g_grow_ages(gpointer data, gpointer user_data) {
    ScheduleSet *x = (ScheduleSet *)data;
    NodeData *    pd = (NodeData *)user_data;

    if (pd->cstat[x->id] == wctlp_LOWER || pd->cstat[x->id] == wctlp_FREE) {
        x->age++;

        if (x->age > pd->retirementage) {
            pd->dzcount++;
        }
    } else {
        x->age = 0;
    }
}

static int grow_ages(NodeData *pd) {
    int val = 0;
    int nb_cols;
    wctlp_get_nb_cols(pd->RMP, &nb_cols);
    assert(nb_cols == pd->localColPool->len);
    CC_IFFREE(pd->cstat, int);
    pd->cstat = (int *)CC_SAFE_MALLOC(nb_cols, int);
    CCcheck_NULL_2(pd->cstat, "Failed to allocate cstat");
    val = wctlp_basis_cols(pd->RMP, pd->cstat, 0);
    CCcheck_val_2(val, "Failed in wctlp_basis_cols");
    pd->dzcount = 0;

    g_ptr_array_foreach(pd->localColPool, g_grow_ages, pd);

CLEAN:
    return val;
}

static int delete_old_cclasses(NodeData *pd) {
    int          val = 0;
    int          min_numdel = pd->njobs * min_ndelrow_ratio;
    int          nb_col;
    guint        i;
    guint        count = pd->localColPool->len;
    ScheduleSet *tmp_schedule;
    /** pd->dzcount can be deprecated! */
    pd->dzcount = 0;

    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet *)g_ptr_array_index(pd->localColPool, i);
        if (tmp_schedule->age > 0) {
            pd->dzcount++;
        }
    }

    if (pd->dzcount > min_numdel) {
        int          it = 0;
        int          first_del = -1;
        int          last_del = -1;
        for (i = 0; i < count; ++i) {
            tmp_schedule = (ScheduleSet *)g_ptr_array_index(pd->localColPool, it);
            if (tmp_schedule->age <= pd->retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = wctlp_deletecols(pd->RMP, first_del, last_del);
                    CCcheck_val_2(val, "Failed in wctlp_deletecols");
                    g_ptr_array_remove_range(pd->localColPool, first_del, last_del - first_del + 1);
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
            g_ptr_array_remove_range(pd->localColPool, first_del, last_del - first_del + 1);
        }

        if (dbg_lvl() > 1) {
            printf("Deleted %d out of %d columns with age > %d.\n", pd->dzcount, count, pd->retirementage);
        }

        wctlp_get_nb_cols(pd->RMP, &nb_col);
        assert(pd->localColPool->len == nb_col);
        for (i = 0; i < pd->localColPool->len; ++i) {
            tmp_schedule = (ScheduleSet *)g_ptr_array_index(pd->localColPool, i);
            tmp_schedule->id = i;
        }
        pd->dzcount = 0;
    }

CLEAN:
    return val;
}

int delete_infeasible_cclasses(NodeData *pd) {
    int          val = 0;
    int          nb_col;
    guint        i;
    guint        count = pd->localColPool->len;
    ScheduleSet *tmp_schedule;
    /** pd->dzcount can be deprecated! */
    pd->dzcount = 0;

    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet *)g_ptr_array_index(pd->localColPool, i);
        if (tmp_schedule->age > 0) {
            pd->dzcount++;
        }
    }

    int          it = 0;
    int          first_del = -1;
    int          last_del = -1;
    for (i = 0; i < count; ++i) {
        tmp_schedule = (ScheduleSet *)g_ptr_array_index(pd->localColPool, it);
        if (tmp_schedule->del != 1) {
            if (first_del != -1) {
                /** Delete recently found deletion range.*/
                val = wctlp_deletecols(pd->RMP, first_del, last_del);
                CCcheck_val_2(val, "Failed in wctlp_deletecols");
                g_ptr_array_remove_range(pd->localColPool, first_del, last_del - first_del + 1);
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
        g_ptr_array_remove_range(pd->localColPool, first_del, last_del - first_del + 1);
    }

    if (dbg_lvl() > 1) {
        printf("Deleted %d out of %d columns with age > %d.\n", pd->dzcount, count, pd->retirementage);
    }

    wctlp_get_nb_cols(pd->RMP, &nb_col);
    assert(pd->localColPool->len == nb_col);
    printf("number of cols = %d\n", nb_col);
    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet *)g_ptr_array_index(pd->localColPool, i);
        tmp_schedule->id = i;
    }
    pd->dzcount = 0;

CLEAN:
    return val;
}

void g_make_pi_feasible(gpointer data, gpointer user_data) {
    ScheduleSet *x = (ScheduleSet *)data;
    NodeData *    pd = (NodeData *)user_data;
    Job *        tmp_j;

    int    i;
    double colsum = .0;

    for (i = 0; i < x->job_list->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(x->job_list, i);
        if (signbit(pd->pi[tmp_j->job])) {
            pd->pi[tmp_j->job] = 0.0;
        }

        colsum += pd->pi[tmp_j->job];
        colsum = nextafter(colsum, DBL_MAX);
    }

    if (!signbit(pd->pi[pd->njobs])) {
        pd->pi[pd->njobs] = 0.0;
    }

    colsum += pd->pi[pd->njobs];
    colsum = nextafter(colsum, DBL_MAX);

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (i = 0; i < x->job_list->len; ++i) {
            tmp_j = (Job *)g_ptr_array_index(x->job_list, i);
            pd->pi[tmp_j->job] /= colsum;
            pd->pi[tmp_j->job] *= x->total_weighted_completion_time;
            newcolsum += pd->pi[tmp_j->job];
        }

        pd->pi[pd->njobs] /= colsum;
        pd->pi[pd->njobs] *= x->total_weighted_completion_time;
        newcolsum += pd->pi[pd->njobs];

        if (dbg_lvl() > 1) {
            printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                   x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible(NodeData *pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible, pd);
}

void g_make_pi_feasible_farkas(gpointer data, gpointer user_data) {
    ScheduleSet *x = (ScheduleSet *)data;
    NodeData *    pd = (NodeData *)user_data;
    Job *        tmp_j;

    int    i;
    double colsum = .0;

    for (i = 0; i < x->job_list->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(x->job_list, i);
        if (signbit(pd->pi[tmp_j->job])) {
            pd->pi[tmp_j->job] = 0.0;
        }

        colsum += pd->pi[tmp_j->job];
        colsum = nextafter(colsum, DBL_MAX);
    }

    colsum += pd->pi[pd->njobs];

    if (colsum > x->total_weighted_completion_time) {
        double newcolsum = .0;
        for (i = 0; i < x->job_list->len; ++i) {
            tmp_j = (Job *)g_ptr_array_index(x->job_list, i);
            pd->pi[tmp_j->job] /= colsum;
            pd->pi[tmp_j->job] *= x->total_weighted_completion_time;
            newcolsum += pd->pi[tmp_j->job];
        }

        pd->pi[pd->njobs] /= colsum;
        pd->pi[pd->njobs] *= x->total_weighted_completion_time;
        newcolsum += pd->pi[pd->njobs];

        if (dbg_lvl() > 1) {
            printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                   x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible_farkas_pricing(NodeData *pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible_farkas, pd);
}

int compute_objective(NodeData *pd, Parms *parms) {
    int val = 0;
    int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    for (i = 0; i < pd->njobs + 1; i++) {
        pd->LP_lower_bound_dual += (double)pd->pi[i] * pd->rhs[i];
    }
    pd->LP_lower_bound_dual -= 0.0001;

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(pd->RMP, &(pd->LP_lower_bound));
    pd->LP_lower_bound -= 0.0001;
    CCcheck_val_2(val, "wctlp_objval failed");
    pd->lower_bound = ((int)ceil(pd->LP_lower_bound_dual) < (int)ceil(pd->LP_lower_bound))
            ? (int)ceil(pd->LP_lower_bound_dual)
            : (int)ceil(pd->LP_lower_bound);
    pd->LP_lower_bound_BB = CC_MIN(pd->LP_lower_bound, pd->LP_lower_bound_dual);
    pd->eta_out = pd->LP_lower_bound_BB;

    if (pd->iterations%pd->njobs == 0) {
        printf(
            "Current primal LP objective: %19.16f  (LP_dual-bound %19.16f, "
            "lowerbound = %d, eta_in = %f, eta_out = %f).\n",
            pd->LP_lower_bound, pd->LP_lower_bound_dual, pd->lower_bound, pd->eta_in, pd->eta_out);
    }

CLEAN:
    return val;
}

int compute_lower_bound(Problem *problem, NodeData *pd) {
    int       j, val = 0;
    int       break_while_loop = 1;
    int       nnonimprovements = 0;
    int       status = GRB_LOADED;
    double    real_time_solve_lp;
    double    real_time_pricing;
    Parms *parms = &(problem->parms);

    // if (pd->status == infeasible) {
    //     wctlp_write(pd->RMP, "test.lp");
    //     wctlp_compute_IIS(pd->RMP);
    //     pd->test = 0;
    //     goto CLEAN;
    // }

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
        add_solution_to_colpool(problem->opt_sol, pd);
    }
    reset_nblayers(pd->jobarray);

    if (!pd->RMP) {
        val = build_rmp(pd, 0);
        CCcheck_val(val, "build_lp failed");
    }

    pd->retirementage = (int)sqrt(pd->njobs) + 30;

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

    /** Compute LP relaxation */
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
            val = wctlp_pi(pd->RMP, pd->pi);
            CCcheck_val_2(val, "wctlp_pi failed");
            /** Compute the objective function */
            val = compute_objective(pd, parms);
            CCcheck_val_2(val, "Failed in compute_objective");
            memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
            pd->eta_out = pd->LP_lower_bound_dual;
            break;

        case WCTLP_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = wctlp_pi_inf(pd->RMP, pd->pi);
            CCcheck_val_2(val, "Failed at wctlp_pi_inf");
            break;
    }

    break_while_loop = 0;
    CCutil_suspend_timer(&(problem->tot_cputime));
    CCutil_resume_timer(&(problem->tot_cputime));

    while ((pd->iterations < pd->maxiterations) && !break_while_loop && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        /** 
         * Delete old columns
         */
        if (pd->dzcount > pd->njobs * min_ndelrow_ratio && status == GRB_OPTIMAL) {
            // val = delete_old_cclasses(pd);
            CCcheck_val_2(val, "Failed in delete_old_cclasses");
        }

        /** 
         * Solve the pricing problem
         */
        real_time_pricing = getRealTime();
        CCutil_start_resume_time(&problem->tot_pricing);

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

                break;
        }

        CCutil_suspend_timer(&problem->tot_pricing);
        real_time_pricing = getRealTime() - real_time_pricing;
        problem->real_time_pricing += real_time_pricing;

        if (pd->update) {
            for (j = 0; j < pd->nnewsets; j++) {
                val = add_scheduleset_to_rmp(pd->newsets + j, pd);
                CCcheck_val_2(val, "wctlp_addcol failed");
                g_ptr_array_add(pd->localColPool, pd->newsets + j);
            }
            pd->newsets = NULL;
        }

        switch (status) {
            case GRB_OPTIMAL:
                switch (parms->stab_technique) {
                    case stab_wentgnes:
                    case stab_dynamic:
                    case stab_hybrid:
                        break_while_loop =
                            (CC_ABS(pd->eta_out - pd->eta_in) < 0.00001) || (ceil(pd->eta_in - 0.00001) >= pd->eta_out);
                        break;

                    case no_stab:
                        break_while_loop =
                            (pd->nnewsets == 0 || nnonimprovements > 5);
                        break;
                }

                break;

            case GRB_INFEASIBLE:
                break_while_loop = (pd->nnewsets == 0);
                break;
        }

        /** Compute LP relaxation */
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
            case GRB_OPTIMAL:
                /** 
                 * grow ages of the different columns
                 */
                val = grow_ages(pd);
                CCcheck_val_2(val, "Failed in grow_ages");

                /** get the dual variables and make them feasible */
                /** Compute the objective function */

                if (pd->update) {
                    val = wctlp_pi(pd->RMP, pd->pi);
                    CCcheck_val_2(val, "wctlp_pi failed");
                    val = compute_objective(pd, parms);
                    CCcheck_val_2(val, "Failed in compute_objective");
                    memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
                }

                break;

            case GRB_INFEASIBLE:
                /** get the dual variables and make them feasible */
                val = wctlp_pi_inf(pd->RMP, pd->pi);
                CCcheck_val_2(val, "wctlp_pi failed");
                break;
        }

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

                if(parms->pricing_solver < dp_solver) {
                    val = wctlp_optimize(pd->RMP, &status);
                    CCcheck_val_2(val, "wctlp_optimize failed");
                    val = compute_objective(pd, parms);
                    CCcheck_val_2(val, "Failed in computing the objective");

                    // reset_nblayers(pd->jobarray);
                    // calculate_nblayers(pd, 2);
                    reduce_cost_fixing(pd);
                    // val = check_schedules(pd);
                    // CCcheck_val_2(val, "Failed in checkschedules");
                    // delete_infeasible_cclasses(pd);
                }
                // pd->status = LP_bound_computed;
                // val = wctlp_pi(pd->RMP, pd->pi);
                // CCcheck_val_2(val, "wctlp_pi failed");

                /**
                 * Compute the objective function
                 */
                val = wctlp_optimize(pd->RMP, &status);
                CCcheck_val_2(val, "wctlp_optimize failed");
                val = compute_objective(pd, parms);
                CCcheck_val_2(val, "Failed in compute_objective");
                memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
                printf("size evolution %lu\n", get_size_graph(pd->solver));
                break;

            case GRB_INFEASIBLE:
                pd->status = infeasible;
                pd->test = 0;
                wctlp_write(pd->RMP, "test1.lp");
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
    problem->global_lower_bound = CC_MAX(pd->lower_bound + pd->problem->off, problem->global_lower_bound);

    if(pd->depth == 0) {
        problem->root_lower_bound = problem->global_lower_bound;
        problem->root_upper_bound = problem->global_upper_bound;
        problem->root_rel_error = (double)(problem->global_upper_bound - problem->global_lower_bound) /
                    ((double)problem->global_lower_bound);
    }

    fflush(stdout);
    problem->nb_generated_col += pd->iterations;
    CCutil_suspend_timer(&(problem->tot_lb));

CLEAN:
    return val;
}

int print_x(NodeData *pd){
    int val = 0;
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status");

    switch(status){
        case GRB_OPTIMAL:
            val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
            CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
            val = wctlp_x(pd->RMP, pd->x, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");

            for(int i = 0; i < nb_cols; ++i) {
                if(pd->x[i] > 0.00001) {
                    printf("x = %f\n", pd->x[i]);
                    ScheduleSet *tmp = (ScheduleSet *) g_ptr_array_index(pd->localColPool,i);
                    int C = 0;
                    for (int j = 0; j < tmp->job_list->len; ++j) {
                        Job *j1 = (Job *) g_ptr_array_index(tmp->job_list, j);
                        
                        if(j  < tmp->job_list->len - 1) {
                            Job *j2 = (Job *) g_ptr_array_index(tmp->job_list, j + 1);
                            if(! bool_diff_Fij(C, j2, j1)) {
                                printf("%d ", bool_diff_Fij(C, j2, j1));
                            }
                        }
                        C += j1->processing_time;
                    }
                    printf("\n");
                    g_ptr_array_foreach(tmp->job_list, g_print_machine, NULL);
                    printf("\n");
                }
            }
        break;
    }

    CLEAN:

    return val;
}

int calculate_nblayers(NodeData *pd, int k){
    int val = 0;
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status");
    if(status == GRB_LOADED) {
        wctlp_optimize(pd->RMP, &status);
    }
    
    reset_nblayers(pd->jobarray);
    
    switch(status){
        case GRB_OPTIMAL:
            val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
            CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
            val = wctlp_x(pd->RMP, pd->x, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");

            for(unsigned i = 0; i < nb_cols; ++i) {
                if(pd->x[i] > 0.00001) {
                    ScheduleSet *tmp = (ScheduleSet *) g_ptr_array_index(pd->localColPool,i);
                        for(int j = 0; j < (int) tmp->job_list->len - k; ++j) {
                            Job* j1 = (Job*) g_ptr_array_index(tmp->job_list, j);
                            Job* j2 = (Job*) g_ptr_array_index(tmp->job_list, j + k);
                            if(j1 == j2) {
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

int calculate_x_e(NodeData *pd){
    int val = 0;
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status")
    
    switch(status){
        case GRB_OPTIMAL:
            val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
            CCcheck_val_2(val, "Failed to get nb cols");
            pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
            CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
            val = wctlp_x(pd->RMP, pd->x, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");
            pd->x_e = CC_SAFE_REALLOC(pd->x_e, 2*get_size_graph(pd->solver), double);
            CCcheck_NULL_2(pd->x_e, "Failed to reallocate memory to  pd->x_e");

            for(unsigned i = 0; i < 2*get_size_graph(pd->solver); ++i) {
                pd->x_e[i] = 0.0;
            }
            construct_lp_sol_from_rmp(pd);
        break;
    }

    CLEAN:

    return val;
}

int check_schedules(NodeData *pd) {
    int val = 0;
    int nb_cols;
    int status;

    val = wctlp_status(pd->RMP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status")
    
    val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
    CCcheck_val_2(val, "Failed to get nb cols");
    assert(nb_cols == pd->localColPool->len);
    printf("number of cols check %d\n", nb_cols);
    for(unsigned i = 0; i < nb_cols; ++i) {
        ScheduleSet *tmp = (ScheduleSet *) g_ptr_array_index(pd->localColPool,i);
        if(check_schedule_set(tmp, pd) == 1) {
            tmp->del = 0;
        } else {
            tmp->del = 1;
            // printf("deleted column\n");
            // print_schedule(tmp, 1);
        }
    }

    if(status != GRB_OPTIMAL) {
        wctlp_compute_IIS(pd->RMP);
        wctlp_write(pd->RMP, "rmpcheck.lp");
        // printf("check file %d\n", status);
        // getchar();
    }

    CLEAN:

    return val;  
}
