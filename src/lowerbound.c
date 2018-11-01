#include <wct.h>
#include <solver.h>

static const double min_ndelrow_ratio = 0.9;

void g_print_ages_col(gpointer data, gpointer user_data) {
    scheduleset *x = (scheduleset *)data;

    printf(" %4d", x->age);
}

/** Help function for column generation */
static void print_ages(wctdata *pd) {
    printf("AGES:");

    g_ptr_array_foreach(pd->localColPool, g_print_ages_col, NULL);

    printf("\n");
}

void g_grow_ages(gpointer data, gpointer user_data) {
    scheduleset *x = (scheduleset *)data;
    wctdata *    pd = (wctdata *)user_data;

    if (pd->cstat[x->id] == wctlp_LOWER || pd->cstat[x->id] == wctlp_FREE) {
        x->age++;

        if (x->age > pd->retirementage) {
            pd->dzcount++;
        }
    } else {
        x->age = 0;
    }
}

static int grow_ages(wctdata *pd) {
    int val = 0;
    int nb_cols;
    wctlp_get_nb_cols(pd->LP, &nb_cols);
    assert(nb_cols == pd->localColPool->len);
    CC_IFFREE(pd->cstat, int);
    pd->cstat = (int *)CC_SAFE_MALLOC(nb_cols, int);
    CCcheck_NULL_2(pd->cstat, "Failed to allocate cstat");
    val = wctlp_basis_cols(pd->LP, pd->cstat, 0);
    CCcheck_val_2(val, "Failed in wctlp_basis_cols");
    pd->dzcount = 0;

    g_ptr_array_foreach(pd->localColPool, g_grow_ages, pd);

CLEAN:
    return val;
}

static int delete_old_cclasses(wctdata *pd) {
    int          val = 0;
    int          it = 0;
    int          min_numdel = pd->njobs * min_ndelrow_ratio;
    int          first_del = -1;
    int          last_del = -1;
    int          nb_col;
    guint        i;
    guint        count = pd->localColPool->len;
    scheduleset *tmp_schedule;
    /** pd->dzcount can be deprecated! */
    pd->dzcount = 0;

    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (scheduleset *)g_ptr_array_index(pd->localColPool, i);
        if (tmp_schedule->age > 0) {
            pd->dzcount++;
        }
    }

    if (pd->dzcount > min_numdel) {
        for (i = 0; i < count; ++i) {
            tmp_schedule = (scheduleset *)g_ptr_array_index(pd->localColPool, it);
            if (tmp_schedule->age <= pd->retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = wctlp_deletecols(pd->LP, first_del, last_del);
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
            wctlp_deletecols(pd->LP, first_del, last_del);
            CCcheck_val_2(val, "Failed in wctlp_deletecols");
            g_ptr_array_remove_range(pd->localColPool, first_del, last_del - first_del + 1);
        }

        if (dbg_lvl() > 1) {
            printf("Deleted %d out of %d columns with age > %d.\n", pd->dzcount, count, pd->retirementage);
        }

        wctlp_get_nb_cols(pd->LP, &nb_col);
        assert(pd->localColPool->len == nb_col);
        for (i = 0; i < pd->localColPool->len; ++i) {
            tmp_schedule = (scheduleset *)g_ptr_array_index(pd->localColPool, i);
            tmp_schedule->id = i;
        }
        pd->dzcount = 0;
    }

CLEAN:
    return val;
}

void g_make_pi_feasible(gpointer data, gpointer user_data) {
    scheduleset *x = (scheduleset *)data;
    wctdata *    pd = (wctdata *)user_data;
    Job *        tmp_j;

    int    i;
    double colsum = .0;

    for (i = 0; i < x->jobs->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(x->jobs, i);
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

    if (colsum > x->totwct) {
        double newcolsum = .0;
        for (i = 0; i < x->jobs->len; ++i) {
            tmp_j = (Job *)g_ptr_array_index(x->jobs, i);
            pd->pi[tmp_j->job] /= colsum;
            pd->pi[tmp_j->job] *= x->totwct;
            newcolsum += pd->pi[tmp_j->job];
        }

        pd->pi[pd->njobs] /= colsum;
        pd->pi[pd->njobs] *= x->totwct;
        newcolsum += pd->pi[pd->njobs];

        if (dbg_lvl() > 1) {
            printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                   x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible(wctdata *pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible, pd);
}

void g_make_pi_feasible_farkas(gpointer data, gpointer user_data) {
    scheduleset *x = (scheduleset *)data;
    wctdata *    pd = (wctdata *)user_data;
    Job *        tmp_j;

    int    i;
    double colsum = .0;

    for (i = 0; i < x->jobs->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(x->jobs, i);
        if (signbit(pd->pi[tmp_j->job])) {
            pd->pi[tmp_j->job] = 0.0;
        }

        colsum += pd->pi[tmp_j->job];
        colsum = nextafter(colsum, DBL_MAX);
    }

    colsum += pd->pi[pd->njobs];

    if (colsum > x->totwct) {
        double newcolsum = .0;
        for (i = 0; i < x->jobs->len; ++i) {
            tmp_j = (Job *)g_ptr_array_index(x->jobs, i);
            pd->pi[tmp_j->job] /= colsum;
            pd->pi[tmp_j->job] *= x->totwct;
            newcolsum += pd->pi[tmp_j->job];
        }

        pd->pi[pd->njobs] /= colsum;
        pd->pi[pd->njobs] *= x->totwct;
        newcolsum += pd->pi[pd->njobs];

        if (dbg_lvl() > 1) {
            printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                   x->id, colsum, newcolsum);
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible_farkas_pricing(wctdata *pd) {
    g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible_farkas, pd);
}

int compute_objective(wctdata *pd, wctparms *parms) {
    int val = 0;
    int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    for (i = 0; i < pd->njobs + 1; i++) {
        pd->LP_lower_bound_dual += (double)pd->pi[i] * pd->rhs[i];
    }
    pd->LP_lower_bound_dual -= 0.0001;

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(pd->LP, &(pd->LP_lower_bound));
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

int compute_lower_bound(wctproblem *problem, wctdata *pd) {
    int       j, val = 0;
    int       break_while_loop = 1;
    int       nnonimprovements = 0;
    int       status = GRB_LOADED;
    int       nb_cols;
    double    real_time_solve_lp;
    double    real_time_pricing;
    wctparms *parms = &(problem->parms);

    if (pd->status == infeasible) {
        goto CLEAN;
    }

    if (dbg_lvl() > 1) {
        printf(
            "Starting compute_lower_bound with lb %d and ub %d at depth %d(id "
            "= "
            "%d, opt_track = %d)\n",
            pd->lower_bound, pd->upper_bound, pd->depth, pd->id, pd->opt_track);
    }

    CCutil_start_resume_time(&(problem->tot_lb));

    /** Construct solutions if list of columns is empty */
    // if (!pd->ccount && parms->construct) {
    //     solution *new_sol;
    //     new_sol = solution_alloc(pd->nmachines, pd->njobs, 0);
    //     construct_edd(problem, new_sol);
    //     partlist_to_scheduleset(new_sol->part, pd->nmachines, pd->njobs,
    //                             &(pd->cclasses), &(pd->ccount));
    //     solution_free(&new_sol);
    //     assert(pd->gallocated >= pd->ccount);
    // }

    if (pd->localColPool->len == 0) {
        solution *new_sol;
        new_sol = solution_alloc(pd->nmachines, pd->njobs, 0);
        construct_edd(problem, new_sol);
        solution_canonical_order(new_sol, pd->local_intervals);
        add_solution_to_colpool(new_sol, pd);
    }
    reset_nblayers(pd->jobarray);

    if (!pd->LP) {
        val = build_lp(pd, 0);
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
    val = wctlp_optimize(pd->LP, &status);
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
            val = wctlp_pi(pd->LP, pd->pi);
            CCcheck_val_2(val, "wctlp_pi failed");
            /** Compute the objective function */
            val = compute_objective(pd, parms);
            CCcheck_val_2(val, "Failed in compute_objective");
            memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
            pd->eta_out = pd->LP_lower_bound_dual;
            break;

        case WCTLP_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = wctlp_pi_inf(pd->LP, pd->pi);
            CCcheck_val_2(val, "Failed at wctlp_pi_inf");
            break;
    }

    break_while_loop = 0;
    CCutil_suspend_timer(&(problem->tot_cputime));
    CCutil_resume_timer(&(problem->tot_cputime));

    while ((pd->iterations < pd->maxiterations) && !break_while_loop && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        /** delete old columns */
        if (pd->dzcount > pd->njobs * min_ndelrow_ratio && status == GRB_OPTIMAL) {
            val = delete_old_cclasses(pd);
        }

        /** Solve the pricing problem*/
        real_time_pricing = getRealTime();
        CCutil_start_resume_time(&problem->tot_pricing);

        switch (status) {
            case GRB_OPTIMAL:
                pd->iterations++;
                pd->status = infeasible;

                if (pd->iterations < pd->maxiterations) {
                    switch (parms->stab_technique) {
                        case stab_wentgnes:
                            if(pd->iterations%5 == 0) {
                                val = solve_pricing(pd, parms,0);
                                CCcheck_val_2(val, "Failed in solving pricing");
                            } else {
                                val = solve_stab(pd, parms);
                                CCcheck_val_2(val, "Failed in solve_stab");
                            }
                            break;

                        case stab_dynamic:
                            val = solve_stab_dynamic(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
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
                val = addColToLP(pd->newsets + j, pd);
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
                        break_while_loop =
                            (CC_OURABS(pd->eta_out - pd->eta_in) < 0.0001) || (ceil(pd->eta_in - 0.0001) >= pd->eta_out);
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
        val = wctlp_optimize(pd->LP, &status);
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
                /** grow ages of the different columns */
                val = grow_ages(pd);
                CCcheck_val_2(val, "Failed in grow_ages");

                /** get the dual variables and make them feasible */
                /** Compute the objective function */

                if (pd->update) {
                    val = wctlp_pi(pd->LP, pd->pi);
                    CCcheck_val_2(val, "wctlp_pi failed");
                    val = compute_objective(pd, parms);
                    CCcheck_val_2(val, "Failed in compute_objective");
                    memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
                }

                break;

            case GRB_INFEASIBLE:
                /** get the dual variables and make them feasible */
                val = wctlp_pi_inf(pd->LP, pd->pi);
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

                /** change status of problem */
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

                wctlp_get_nb_cols(pd->LP, &nb_cols);
                pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
                CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
                val = wctlp_x(pd->LP, pd->x, 0);
                CCcheck_val_2(val, "Failed in wctlp_x");
                calculate_nblayers(pd);
                val = calculate_x_e(pd);
                CCcheck_val_2(val, "Failed in calculate_x_e");
                pd->status = LP_bound_computed;
                val = wctlp_pi(pd->LP, pd->pi);
                CCcheck_val_2(val, "wctlp_pi failed");
                solve_pricing(pd, parms, 1);
                if(pd->nnewsets) {
                    schedulesets_free(&pd->newsets, &(pd->nnewsets));
                }
                /** Compute the objective function */
                val = compute_objective(pd, parms);
                CCcheck_val_2(val, "Failed in compute_objective");
                memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
                evaluate_nodes(pd);
                calculate_new_ordered_jobs(pd);
                break;

            case GRB_INFEASIBLE:
                pd->status = infeasible;
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
    }
    problem->global_lower_bound = CC_MAX(pd->lower_bound + pd->problem->off, problem->global_lower_bound);

    fflush(stdout);
    problem->nb_generated_col += pd->iterations;
    CCutil_suspend_timer(&(problem->tot_lb));

CLEAN:
    return val;
}

int print_x(wctdata *pd){
    int val = 0;
    int nb_cols;
    int status;

    wctlp_get_nb_cols(pd->LP, &nb_cols);
    pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    val = wctlp_x(pd->LP, pd->x, 0);
    CCcheck_val_2(val, "Failed in wctlp_x");
    wctlp_status(pd->LP, &status);
    switch(status){
        case GRB_OPTIMAL:
            for(unsigned i = 0; i < nb_cols; ++i) {
                if(pd->x[i] > 0.0001) {
                    printf("x = %f\n", pd->x[i]);
                    scheduleset *tmp = (scheduleset *) g_ptr_array_index(pd->localColPool,i);
                    g_ptr_array_foreach(tmp->jobs, g_print_machine, NULL);
                    printf("\n");
                }
            }
        break;
    }
    CLEAN:
    return val;
}

int calculate_nblayers(wctdata *pd){
    int val = 0;
    int nb_cols;
    int status;

    wctlp_get_nb_cols(pd->LP, &nb_cols);
    pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
    reset_nblayers(pd->jobarray);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    val = wctlp_x(pd->LP, pd->x, 0);
    CCcheck_val_2(val, "Failed in wctlp_x");
    wctlp_status(pd->LP, &status);
    switch(status){
        case GRB_OPTIMAL:
            for(unsigned i = 0; i < nb_cols; ++i) {
                if(pd->x[i] > 0.0001) {
                    scheduleset *tmp = (scheduleset *) g_ptr_array_index(pd->localColPool,i);
                    g_ptr_array_foreach(tmp->jobs, g_compute_nblayers_schedule, tmp);
                }
            }
        break;
    }
    CLEAN:
    return val;
}

int calculate_x_e(wctdata *pd){
    int val = 0;
    int nb_cols;
    int status;

    wctlp_get_nb_cols(pd->LP, &nb_cols);
    pd->x = CC_SAFE_REALLOC(pd->x, nb_cols, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    val = wctlp_x(pd->LP, pd->x, 0);
    CCcheck_val_2(val, "Failed in wctlp_x");
    wctlp_status(pd->LP, &status);
    switch(status){
        case GRB_OPTIMAL:
            for(unsigned i = 0; i < 2*get_datasize(pd->solver); ++i) {
                pd->x_e[i] = 0.0;
            }
            for(unsigned i = 0; i < nb_cols; ++i) {
                if(pd->x[i] > 0.00001) {
                    scheduleset *tmp = (scheduleset *) g_ptr_array_index(pd->localColPool,i);
                    for(unsigned j = 0; j < tmp->e_list->len; ++j) {
                        int *ptr = (int *) g_ptr_array_index(tmp->e_list, j);
                        pd->x_e[*ptr] += pd->x[i];
                    }
                }
            }
        break;
    }

    CLEAN:
    return val;
}
