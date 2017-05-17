#include <wct.h>

static const double min_ndelrow_ratio = 0.3;

/** Help function for column generation */
static void print_ages(wctdata *pd) {
    int i;
    printf("AGES:");

    for (i = 0; i < pd->ccount; ++i) {
        printf(" %4d", pd->cclasses[i].age);
    }

    printf("\n");
}

static int grow_ages(wctdata *pd) {
    int  val = 0;
    int  i;
    int *cstat;
    cstat = (int *)CC_SAFE_MALLOC(pd->ccount, int);
    CCcheck_NULL_2(cstat, "Failed to allocate cstat");
    val = wctlp_basis_cols(pd->LP, cstat, 0);
    CCcheck_val_2(val, "Failed in wctlp_basis_cols");
    pd->dzcount = 0;

    for (i = 0; i < pd->ccount; ++i) {
        if (cstat[i] == wctlp_LOWER || cstat[i] == wctlp_FREE) {
            pd->cclasses[i].age++;

            if (pd->cclasses[i].age > pd->retirementage) {
                pd->dzcount++;
            }
        } else {
            pd->cclasses[i].age = 0;
        }
    }

/*    printf("%d out of %d are older than %d.\n", pd->dzcount, pd->ccount,  */
/*           pd->retirementage); */
CLEAN:
    CC_IFFREE(cstat, int);
    return val;
}

MAYBE_UNUSED
static int delete_old_cclasses(wctdata *pd) {
    int val = 0;
    int i;
    int min_numdel = pd->njobs * min_ndelrow_ratio;
    int first_del = -1;
    int last_del = -1;
    /** pd->dzcount can be deprecated! */
    pd->dzcount = 0;

    for (i = 0; i < pd->ccount; ++i) {
        if (pd->cclasses[i].age > pd->retirementage) {
            pd->dzcount++;
        }
    }

    if (pd->dzcount > min_numdel) {
        int          new_ccount = 0;
        scheduleset *new_cclasses = CC_SAFE_MALLOC(pd->gallocated, scheduleset);
        CCcheck_NULL_2(new_cclasses, "Failed to allocate new_cclasses");
        assert(pd->gallocated >= pd->ccount);

        for (i = 0; i < pd->gallocated; ++i) {
            scheduleset_init(new_cclasses + i);
        }

        for (i = 0; i < pd->ccount; ++i) {
            if (pd->cclasses[i].age <= pd->retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = wctlp_deletecols(pd->LP, first_del, last_del);
                    CCcheck_val_2(val, "Failed in wctlp_deletecols");
                    first_del = last_del = -1;
                }

                memcpy(new_cclasses + new_ccount, pd->cclasses + i,
                       sizeof(scheduleset));
                new_ccount++;
            } else {
                scheduleset_free(pd->cclasses + i);

                if (first_del == -1) {
                    first_del = new_ccount;
                    last_del = first_del;
                } else {
                    last_del++;
                }
            }
        }

        if (first_del != -1) {
            /** Delete the final range. This can occur if the last
             element is to be deleted, e.g. when no further columns were
             added in a B&B branch.
             */
            wctlp_deletecols(pd->LP, first_del, last_del);
            CCcheck_val_2(val, "Failed in wctlp_deletecols");
        }

        assert(pd->dzcount == pd->ccount - new_ccount);
        CC_IFFREE(pd->cclasses, scheduleset);
        pd->cclasses = new_cclasses;
        pd->ccount = new_ccount;

        if (dbg_lvl() > 1) {
            printf("Deleted %d out of %d columns with age > %d.\n", pd->dzcount,
                   pd->dzcount + pd->ccount, pd->retirementage);
        }

        pd->dzcount = 0;
    }

CLEAN:
    return val;
}

static void reset_ages(scheduleset *cclasses, int cccount) {
    int i;

    for (i = 0; i < cccount; i++) {
        cclasses[i].age = 0;
    }
}

int add_newsets(wctdata *pd) {
    int          val = 0;
    scheduleset *tmpsets = (scheduleset *)NULL;
    int          i;

    if (pd->nnewsets == 0) {
        return val;
    }

    reset_ages(pd->newsets, pd->nnewsets);

    if (pd->ccount + pd->nnewsets > pd->gallocated) {
        pd->gallocated *= 2;
        tmpsets = CC_SAFE_MALLOC(pd->gallocated, scheduleset);
        CCcheck_NULL_2(tmpsets, "Failed to allocate memory to tmpsets");
        memcpy(tmpsets, pd->cclasses, pd->ccount * sizeof(scheduleset));
        free(pd->cclasses);
        pd->cclasses = tmpsets;
        tmpsets = NULL;
    }

    memcpy(pd->cclasses + pd->ccount, pd->newsets,
           pd->nnewsets * sizeof(scheduleset));
    pd->ccount += pd->nnewsets;

    for (i = pd->ccount; i < pd->gallocated; i++) {
        scheduleset_init(pd->cclasses + i);
    }

CLEAN:

    if (val) {
        CC_IFFREE(pd->cclasses, scheduleset);
    }

    CC_IFFREE(pd->newsets, scheduleset);
    pd->nnewsets = 0;
    return val;
}

MAYBE_UNUSED
void make_pi_feasible(wctdata *pd) {
    int c;

    for (c = 0; c < pd->ccount; ++c) {
        int    i;
        double colsum = .0;

        for (i = 0; i < pd->cclasses[c].count; ++i) {
            if (signbit(pd->pi[pd->cclasses[c].members[i]])) {
                pd->pi[pd->cclasses[c].members[i]] = 0.0;
            }

            colsum += pd->pi[pd->cclasses[c].members[i]];
            colsum = nextafter(colsum, DBL_MAX);
        }

        if (!signbit(pd->pi[pd->cclasses[c].members[i]])) {
            pd->pi[pd->cclasses[c].members[i]] = 0.0;
        }

        colsum += pd->pi[pd->cclasses[c].members[i]];
        colsum = nextafter(colsum, DBL_MAX);

        if (colsum > pd->cclasses[c].totwct) {
            double newcolsum = .0;
            for (i = 0; i < pd->cclasses[c].count; ++i) {
                pd->pi[pd->cclasses[c].members[i]] /= colsum;
                pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
                newcolsum += pd->pi[pd->cclasses[c].members[i]];
            }

            pd->pi[pd->cclasses[c].members[i]] /= colsum;
            pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
            newcolsum += pd->pi[pd->cclasses[c].members[i]];

            if (dbg_lvl() > 1) {
                printf(
                    "Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                    c, colsum, newcolsum);
            }
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible_farkas_pricing(wctdata *pd) {
    int c;

    for (c = 0; c < pd->ccount; ++c) {
        int    i;
        double colsum = .0;

        for (i = 0; i < pd->cclasses[c].count; ++i) {
            if (signbit(pd->pi[pd->cclasses[c].members[i]])) {
                pd->pi[pd->cclasses[c].members[i]] = 0.0;
            }

            colsum += pd->pi[pd->cclasses[c].members[i]];
            colsum = nextafter(colsum, DBL_MAX);
        }

        colsum += pd->pi[pd->cclasses[c].members[i]];

        if (colsum > pd->cclasses[c].totwct) {
            double newcolsum = .0;
            for (i = 0; i < pd->cclasses[c].count; ++i) {
                pd->pi[pd->cclasses[c].members[i]] /= colsum;
                pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
                newcolsum += pd->pi[pd->cclasses[c].members[i]];
            }

            pd->pi[pd->cclasses[c].members[i]] /= colsum;
            pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
            newcolsum += pd->pi[pd->cclasses[c].members[i]];

            if (dbg_lvl() > 1) {
                printf(
                    "Decreased column sum of %5d from  %30.20f to  %30.20f\n",
                    c, colsum, newcolsum);
            }
        }
    }
}

int compute_objective(wctdata *pd, wctparms *parms) {
    int val = 0;
    int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    for (i = 0; i < pd->njobs + 1; i++) {
        pd->LP_lower_bound_dual += (double)pd->pi[i] * pd->rhs[i];
    }

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(pd->LP, &(pd->LP_lower_bound));
    CCcheck_val_2(val, "wctlp_objval failed");
    pd->lower_bound =
        ((int)ceil(pd->LP_lower_bound_dual) < (int)ceil(pd->LP_lower_bound))
            ? (int)ceil(pd->LP_lower_bound_dual)
            : (int)ceil(pd->LP_lower_bound);
    pd->LP_lower_bound_BB = CC_MIN(pd->LP_lower_bound, pd->LP_lower_bound_dual);

    if (parms->stab_technique == stab_wentgnes ||
        parms->stab_technique == stab_dynamic) {
        pd->lower_bound = (int)ceil(pd->eta_in - 0.001);
        pd->LP_lower_bound_BB = CC_MIN(pd->LP_lower_bound_BB, pd->eta_in);
    }

    if (dbg_lvl() > 0) {
        printf(
            "Current primal LP objective: %19.16f  (LP_dual-bound %19.16f, "
            "lowerbound = %d).\n",
            pd->LP_lower_bound, pd->LP_lower_bound_dual, pd->lower_bound);
    }

CLEAN:
    return val;
}

int compute_lower_bound(wctproblem *problem, wctdata *pd) {
    int       j, val = 0;
    int       iterations = 0;
    int       break_while_loop = 1;
    int       nnonimprovements = 0;
    int       status = GRB_LOADED;
    double    cur_cputime;
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

    CCutil_start_resume_time(&(problem->tot_lb_lp));

    /** Construct solutions if list of columns is empty */
    if (!pd->ccount && parms->construct) {
        solution *new_sol;
        new_sol = solution_alloc(pd->nmachines, pd->njobs, 0);
        construct_edd(problem, new_sol);
        partlist_to_scheduleset(new_sol->part, pd->nmachines, pd->njobs,
                                &(pd->cclasses), &(pd->ccount));
        solution_free(&new_sol);
        assert(pd->gallocated >= pd->ccount);
    }

    if (!pd->LP) {
        val = build_lp(pd, parms->construct);
        CCcheck_val(val, "build_lp failed");
    }

    pd->retirementage = (int)sqrt(pd->njobs) + 30;

    /** Init alpha */
    switch (parms->stab_technique) {
        case stab_wentgnes:
            pd->alpha = 0.8;
            break;

        case stab_dynamic:
            pd->alpha = 0.0;
            break;

        case no_stab:
            break;
    }

    /** Compute LP relaxation */
    cur_cputime = CCutil_zeit();
    val = wctlp_optimize(pd->LP, &status);
    CCcheck_val_2(val, "wctlp_optimize failed");
    cur_cputime = CCutil_zeit() - cur_cputime;

    if (dbg_lvl() > 1) {
        printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
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
            val = wctlp_pi(pd->LP, pd->pi);
            CCcheck_val_2(val, "wctlp_pi failed");
            /** Compute the objective function */
            val = compute_objective(pd, parms);
            CCcheck_val_2(val, "Failed in compute_objective");
            memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
            pd->eta_out = pd->LP_lower_bound_dual;
            break;

        case GRB_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = wctlp_pi_inf(pd->LP, pd->pi);
            CCcheck_val_2(val, "Failed at wctlp_pi_inf");
            break;
    }

    break_while_loop = 0;
    CCutil_suspend_timer(&(problem->tot_cputime));
    CCutil_resume_timer(&(problem->tot_cputime));

    while ((iterations < pd->maxiterations) && !break_while_loop &&
           problem->tot_cputime.cum_zeit <=
               problem->parms.branching_cpu_limit) {
        /** delete old columns */
        if (pd->dzcount > pd->njobs * min_ndelrow_ratio &&
            status == GRB_OPTIMAL) {
            val = delete_old_cclasses(pd);
        }

        /** Solve the pricing problem*/
        CCutil_start_resume_time(&problem->tot_pricing);

        switch (status) {
            case GRB_OPTIMAL:
                iterations++;
                pd->status = infeasible;

                if (iterations < pd->maxiterations) {
                    switch (parms->stab_technique) {
                        case stab_wentgnes:
                            val = solve_stab(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case stab_dynamic:
                            val = solve_stab_dynamic(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case no_stab:
                            val = solve_pricing(pd, parms);
                            CCcheck_val_2(val, "Failed in solving pricing");
                            break;
                    }
                }

                break;

            case GRB_INFEASIBLE:
                switch (parms->solver) {
                    case bdd_solver:
                    case zdd_solver:
                        val = solve_farkas_dbl(pd);
                        CCcheck_val_2(val, "Failed in solving farkas");
                        break;

                    case DP_solver:
                        val = solve_farkas_dbl_DP(pd);
                        CCcheck_val_2(val, "Failed in solving farkas DP");
                }

                break;
        }

        CCutil_suspend_timer(&problem->tot_pricing);

        for (j = 0; j < pd->nnewsets && pd->update; j++) {
            val = wctlp_addcol(pd->LP, pd->newsets[j].count + 1,
                               pd->newsets[j].members, pd->coef,
                               pd->newsets[j].totwct, 0.0, GRB_INFINITY,
                               wctlp_CONT, NULL);
            CCcheck_val_2(val, "wctlp_addcol failed");
        }

        switch (status) {
            case GRB_OPTIMAL:
                switch (parms->stab_technique) {
                    case stab_wentgnes:
                    case stab_dynamic:
                        break_while_loop =
                            (CC_OURABS(pd->eta_out - pd->eta_in) < 0.00001);
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

        if (pd->update) {
            add_newsets(pd);
        } else {
            schedulesets_free(&pd->newsets, &pd->nnewsets);
        }

        /** Compute LP relaxation */
        cur_cputime = CCutil_zeit();
        val = wctlp_optimize(pd->LP, &status);
        CCcheck_val_2(val, "wctlp_optimize failed");
        cur_cputime = CCutil_zeit() - cur_cputime;

        if (dbg_lvl() > 1) {
            printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
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
                    memcpy(pd->pi_out, pd->pi,
                           sizeof(double) * (pd->njobs + 1));
                    pd->eta_out = pd->LP_lower_bound_dual < pd->LP_lower_bound
                                      ? pd->LP_lower_bound
                                      : pd->LP_lower_bound_dual;
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

    if (iterations < pd->maxiterations &&
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
                        pd->id, iterations, pd->opt_track);
                }

                pd->x = CC_SAFE_MALLOC(pd->ccount, double);
                CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
                val = wctlp_x(pd->LP, pd->x, 0);
                CCcheck_val_2(val, "Failed in wctlp_x");
                pd->status = LP_bound_computed;
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

    if (dbg_lvl() > 1) {
        printf("iterations = %d\n", iterations);
    }

    fflush(stdout);
    problem->nb_generated_col += iterations;
    CCutil_suspend_timer(&(problem->tot_lb_lp));
CLEAN:
    return val;
}
