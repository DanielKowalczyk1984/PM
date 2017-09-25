#include <wct.h>

int grab_int_sol(wctdata *pd, double *x, double tolerance)
{
    int    val = 0;
    double test_incumbent = .0;
    double incumbent;
    int    i;
    int    tot_weighted = 0;
    int   *colored = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(colored, "Failed to allocate colored");
    fill_int(colored, pd->njobs, 0);
    val = wctlp_objval(pd->LP, &incumbent);
    CCcheck_val_2(val, "wctlp_objval failed");
    schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    pd->bestcolors = CC_SAFE_MALLOC(pd->nmachines, scheduleset);
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    pd->nbbest = 0;

    for (i = 0; i < pd->ccount; ++i) {
        test_incumbent += x[i];

        if (x[i] >= 1.0 - tolerance) {
            int j = pd->nbbest;
            int k;
            scheduleset_init(pd->bestcolors + j);
            pd->bestcolors[j].members =
                CC_SAFE_MALLOC(pd->cclasses[i].count, int);
            CCcheck_NULL_2(pd->bestcolors[j].members,
                           "Failed to realloc pd->bestcolors[j].members");

            for (k = 0; k < pd->cclasses[i].count; ++k) {
                if (!colored[pd->cclasses[i].members[k]]) {
                    colored[pd->cclasses[i].members[k]] = 1;
                    pd->bestcolors[j].members[pd->bestcolors[j].count++] =
                        pd->cclasses[i].members[k];
                    pd->bestcolors[j].totweight +=
                        pd->jobarray[pd->cclasses[i].members[k]].processingime;
                    pd->bestcolors[j].totwct +=
                        pd->jobarray[pd->cclasses[i].members[k]].weight *
                        pd->bestcolors[j].totweight;
                }
            }

            pd->nbbest++;
            tot_weighted += pd->bestcolors[j].totwct;

            if (pd->nbbest > pd->nmachines) {
                printf(
                    "ERROR: \"Integral\" solution turned out to be not "
                    "integral!\n");
                fflush(stdout);
                val = 1;
                goto CLEAN;
            }
        }
    }

    val = scheduleset_check(pd->bestcolors, pd->nbbest, pd->njobs);
    CCcheck_val_2(val, "ERROR: An incorrect coloring was created.");
    printf("Intermediate schedule:\n");
    print_schedule(pd->bestcolors, pd->nbbest);
    printf("with total weight %d\n", tot_weighted);
    assert(fabs((double)tot_weighted - incumbent) <= 0.00001);

    if (tot_weighted < pd->upper_bound) {
        pd->upper_bound = tot_weighted;
        pd->besttotwct = tot_weighted;
    }

    if (pd->upper_bound == pd->lower_bound) {
        pd->status = finished;
    }

CLEAN:
    CC_IFFREE(colored, int);
    return val;
}

int addColToLP(scheduleset *set, wctdata *pd){
    wctproblem *problem = pd->problem;
    int val = 0;
    int cind;
    int vind;
    double cval;
    int nb_intervals = (int) pd->local_intervals->len;
    int u = set->nb_interval;
    int njobs = problem->njobs;
    GPtrArray *members = set->partial_machine;
    wctlp *lp = pd->LP;
    Job *job;
    interval *I = (interval *) g_ptr_array_index(pd->local_intervals, u);

    val = wctlp_get_nb_cols(lp, &(set->id_right));
    CCcheck_val_2(val, "Failed to get the number of cols");
    vind = set->id_right;
    val = wctlp_addcol(lp, 0, NULL, NULL, (double) set->totwct_right, 0.0, GRB_INFINITY, wctlp_CONT, NULL);
    CCcheck_val_2(val, "Failed to add column to lp")

    for(unsigned i = 0; i < members->len; ++i) {
        job = (Job *) g_ptr_array_index(members, i);
        cval = 1.0;
        wctlp_chgcoef(lp, 1, &(job->job), &(set->id_right), &cval);
    }

    if(u == 0) {
        cind =  njobs + u;
        cval = (double) -set->totweight + set->c_right + I->b - I->a;
        wctlp_chgcoef(lp, 1, &cind, &vind, &cval);
    } else if (u < nb_intervals - 1){
        int tcind[] = {njobs + u, njobs + u - 1};
        int tvind[] = {vind, vind};
        double tcval[] = {(double) -set->totweight + set->c_right + I->b - I->a, -(double) set->c_right};
        wctlp_chgcoef(lp, 2, tcind, tvind, tcval);
    } else {
        cind = njobs + u - 1;
        cval = -(double) set->c_right;
        wctlp_chgcoef(lp, 1, &cind, &vind, &cval);
    }

    cind =  njobs + nb_intervals - 1 + u;
    cval = 1.0;
    wctlp_chgcoef(lp, 1, &cind, &vind, &cval);

    if(u != 0) {
        val = wctlp_get_nb_cols(lp, &(set->id_left));
        CCcheck_val_2(val, "Failed to get the number of cols");
        vind = set->id_left;
        val = wctlp_addcol(lp, 0, NULL, NULL, (double) set->totwct_left, 0.0, GRB_INFINITY, wctlp_CONT, NULL);
        CCcheck_val_2(val, "Failed to add column to lp")

        for(unsigned i = 0; i < members->len; ++i) {
            job = (Job *) g_ptr_array_index(members, i);
            cval = 1.0;
            wctlp_chgcoef(lp, 1, &(job->job), &(set->id_left), &cval);
        }

        if (u < nb_intervals - 1){
            int tcind[] = {njobs + u, njobs + u - 1};
            int tvind[] = {vind, vind};
            double tcval[] = {(double) -set->totweight + set->c_left + I->b - I->a, -(double) set->c_left};
            wctlp_chgcoef(lp, 2, tcind, tvind, tcval);
        } else {
            cind = njobs + u - 1;
            cval = -(double) set->c_left;
            wctlp_chgcoef(lp, 1, &cind, &vind, &cval);
        }

        cind = njobs + nb_intervals - 1 + u;
        cval = 1.0;
        wctlp_chgcoef(lp, 1, &cind, &vind, &cval);
    }



    CLEAN:
    return val;
}

int build_lp(wctdata *pd, int construct)
{
    int  val = 0;
    wctproblem *problem = pd->problem;
    int  njobs = problem->njobs;
    int nmachines = problem->nmachines;
    int  counter = 0;
    int nb_row;
    int *u = CC_SAFE_MALLOC(pd->local_intervals->len, int);
    fill_int(u, pd->local_intervals->len, 0);
    GPtrArray *intervals = pd->local_intervals;
    interval *I;
    int *covered = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(covered, "Failed to allocate memory to covered");
    val = wctlp_init(&(pd->LP), NULL);
    CCcheck_val_2(val, "wctlp_init failed");
    fill_int(covered, njobs, 0);

    for (int i = 0; i < njobs; i++) {
        val = wctlp_addrow(pd->LP, 0, (int *)NULL, (double *)NULL,
                           wctlp_GREATER_EQUAL, 1.0, (char *)NULL);
        CCcheck_val_2(val, "Failed wctlp_addrow");
    }

    for (unsigned i = 0; i < intervals->len - 1; ++i) {
        I = (interval *) g_ptr_array_index(intervals, i);
        val = wctlp_addrow(pd->LP, 0, (int *)NULL, (double *)NULL,
                           wctlp_EQUAL, 0.0, (char *)NULL);
    }

    for (unsigned i = 0; i < intervals->len; ++i) {
        val = wctlp_addrow(pd->LP, 0, (int *)NULL, (double *)NULL,
                           wctlp_EQUAL, (double) nmachines, (char *)NULL);
    }

    wctlp_get_nb_rows(pd->LP, &nb_row);

    for(unsigned i = problem->nArtificials; i < problem->ColPool->len; ++i) {
        scheduleset *tmp = (scheduleset *) g_ptr_array_index(problem->ColPool, i);
        u[tmp->nb_interval] = 1;
    }

    /** add artificial columns */
    for(int j = 0; j < problem->nArtificials; ++j) {
        scheduleset *tmp = (scheduleset *) g_ptr_array_index(problem->ColPool, j);
        addColToLP(tmp, pd);
        g_ptr_array_add(pd->localColPool, tmp);
    }

    for(unsigned i = problem->nArtificials; i < problem->ColPool->len; ++i) {
        scheduleset *tmp = (scheduleset *) g_ptr_array_index(problem->ColPool, i);
        addColToLP(tmp, pd);
        g_ptr_array_add(pd->localColPool, tmp);
    }

    wctlp_optimize(pd->LP, &counter);
    printf("test off %d\n", problem->off);
    wctlp_get_nb_cols(pd->LP, &counter);
    pd->x = CC_SAFE_MALLOC(counter, double);
    wctlp_x(pd->LP, pd->x, 0);

    for(unsigned i = problem->nArtificials; i < pd->localColPool->len; ++i) {
        scheduleset *tmp = (scheduleset *) g_ptr_array_index(pd->localColPool, i);
        printf("interval %d\n",tmp->nb_interval);
        if(tmp->nb_interval == 0) {
            printf("%f w=%d, p= %d, c=%d jobs= ", pd->x[tmp->id_right], tmp->totwct_right, tmp->totweight, tmp->c_right);
            g_ptr_array_foreach(tmp->partial_machine, g_print_job, NULL);
            printf("\n");
        } else {
            printf("%f w=%d, p= %d, c=%d\n", pd->x[tmp->id_right], tmp->totwct_right, tmp->totweight, tmp->c_right);
            printf("%f w=%d, p= %d, c=%d jobs= ", pd->x[tmp->id_left], tmp->totwct_left, tmp->totweight, tmp->c_left);
            g_ptr_array_foreach(tmp->partial_machine, g_print_job, NULL);
            printf("\n");
        }
    }
    wctlp_write(pd->LP, "test.lp");

    /** add constraint about number of machines */
    // if (construct) {
    //     for (i = 0; i < pd->ccount; i++) {
    //         val = wctlp_addcol(pd->LP, pd->cclasses[i].count + 1,
    //                            pd->cclasses[i].members, pd->coef,
    //                            pd->cclasses[i].totwct, 0.0, GRB_INFINITY,
    //                            wctlp_CONT, NULL);
    //         if (val) {
    //             wctlp_printerrorcode(val);
    //         }
    //         for (j = 0; j < pd->cclasses[i].count && counter < njobs; j++) {
    //             if (!covered[pd->cclasses[i].members[j]]) {
    //                 covered[pd->cclasses[i].members[j]] = 1;
    //                 counter++;
    //             }
    //         }
    //     }
    // }
    // if (counter < njobs) {
    //     /** Farkas Pricing */
    //     for (i = 0; i < njobs; i++) {
    //         if (!covered[i]) {
    //             pd->nnewsets = 1;
    //             pd->newsets = CC_SAFE_MALLOC(1, scheduleset);
    //             scheduleset_init(pd->newsets);
    //             pd->newsets[0].members = CC_SAFE_MALLOC(2, int);
    //             CCcheck_NULL_2(
    //                 pd->newsets[0].members,
    //                 "Failed to allocate memory to pd->newsets->members");
    //             pd->newsets[0].C = CC_SAFE_MALLOC(1, int);
    //             CCcheck_NULL_2(pd->newsets[0].C, "Failed to allocate memory");
    //             pd->newsets[0].table =
    //                 g_hash_table_new(g_direct_hash, g_direct_equal);
    //             CCcheck_NULL_2(pd->newsets[0].table,
    //                            "Failed to allocate memory");
    //             pd->newsets[0].count++;
    //             pd->newsets[0].members[0] = i;
    //             pd->newsets[0].members[1] = njobs;
    //             pd->newsets[0].totwct =
    //                 pd->jobarray[i].weight * pd->jobarray[i].processingime;
    //             pd->newsets[0].totweight = pd->jobarray[i].processingime;
    //             pd->newsets[0].C[0] = pd->jobarray[i].processingime;
    //             g_hash_table_insert(pd->newsets[0].table, GINT_TO_POINTER(i),
    //                                 pd->newsets[0].C);
    //             pd->newsets->age = 0;
    //             val = wctlp_addcol(pd->LP, 2, pd->newsets[0].members, pd->coef,
    //                                pd->newsets[0].totwct, 0.0, GRB_INFINITY,
    //                                wctlp_CONT, NULL);
    //             CCcheck_val_2(val, "Failed in wctlp_addcol");
    //             if (pd->gallocated == 0 && pd->ccount == 0) {
    //                 pd->cclasses = CC_SAFE_MALLOC(njobs, scheduleset);
    //                 for (j = 0; j < njobs; ++j) {
    //                     scheduleset_init(pd->cclasses + i);
    //                 }
    //                 pd->gallocated = njobs;
    //             }
    //             add_newsets(pd);
    //         }
    //     }
    // }
    pd->pi = (double *)CCutil_reallocrus(pd->pi, nb_row * sizeof(double));
    CCcheck_NULL_2(pd->pi, "Failed to allocate memory to pd->pi");
    pd->pi_in = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->pi_in, "Failed to allocate memory");
    fill_dbl(pd->pi_in, nb_row, 0.0);
    pd->eta_in = 0.0;
    pd->pi_out = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->pi_out, "Failed to allocate memory");
    pd->eta_out = 0.0;
    fill_dbl(pd->pi_out, nb_row, 0.0);
    pd->pi_sep = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->pi_sep, "Failed to allocate memory");
    pd->subgradient_in = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->subgradient_in, "Failed to allocate memory");
    pd->subgradient = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->subgradient, "Failed to allocate memory");
    pd->rhs = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->rhs, "Failed to allocate memory");
    val = wctlp_get_rhs(pd->LP, pd->rhs);
    CCcheck_val_2(val, "Failed to get RHS");
CLEAN:

    if (val) {
        wctlp_free(&(pd->LP));
        CC_IFFREE(pd->coef, double);
        CC_IFFREE(pd->pi, double);
        CC_IFFREE(pd->pi_in, double)
        CC_IFFREE(pd->pi_out, double)
        CC_IFFREE(pd->pi_sep, double)
        CC_IFFREE(pd->subgradient, double)
        CC_IFFREE(pd->subgradient_in, double)
        CC_IFFREE(pd->rhs, double)
    }

    CC_IFFREE(covered, int);
    CC_IFFREE(u, int);
    return val;
}
