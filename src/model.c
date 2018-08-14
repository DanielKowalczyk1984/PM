#include <wct.h>

int grab_int_sol(wctdata *pd, double *x, double tolerance) {
    int          val = 0;
    double       test_incumbent = .0;
    double       incumbent;
    int          i;
    int          tot_weighted = 0;
    int          nb_cols;
    scheduleset *tmp_schedule;
    Job *        tmp_j;

    val = wctlp_objval(pd->LP, &incumbent);
    CCcheck_val_2(val, "wctlp_objval failed");
    val = wctlp_get_nb_cols(pd->LP, &nb_cols);
    CCcheck_val_2(val, "Failed get nb_cols");

    schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    pd->bestcolors = CC_SAFE_MALLOC(pd->nmachines, scheduleset);
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    pd->nbbest = 0;

    assert(nb_cols == pd->localColPool->len);
    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (scheduleset *)g_ptr_array_index(pd->localColPool, i);
        test_incumbent += x[i];

        if (x[i] >= 1.0 - tolerance) {
            int j = pd->nbbest;
            int k;
            scheduleset_init(pd->bestcolors + j);

            g_ptr_array_set_size(pd->bestcolors[j].jobs,
                                 tmp_schedule->jobs->len);
            for (k = 0; k < tmp_schedule->jobs->len; ++k) {
                tmp_j = (Job *)g_ptr_array_index(tmp_schedule->jobs, k);
                g_ptr_array_add(pd->bestcolors[j].jobs, tmp_j);
                pd->bestcolors[j].totweight += tmp_j->processingime;
                pd->bestcolors[j].totwct +=
                    value_Fj(pd->bestcolors[j].totweight, tmp_j);
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

    /** Write a check function */
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
    return val;
}

int addColToLP(scheduleset *set, wctdata *pd) {
    int        val = 0;
    int        cind;
    int        vind;
    double     cval;
    int        njobs = pd->njobs;
    GPtrArray *members = set->jobs;
    wctlp *    lp = pd->LP;
    Job *      job;

    val = wctlp_get_nb_cols(lp, &(set->id));
    CCcheck_val_2(val, "Failed to get the number of cols");
    vind = set->id;
    val = wctlp_addcol(lp, 0, NULL, NULL, (double)set->totwct, 0.0,
                       GRB_INFINITY, wctlp_CONT, NULL);
    CCcheck_val_2(val, "Failed to add column to lp")

        for (unsigned i = 0; i < members->len; ++i) {
        job = (Job *)g_ptr_array_index(members, i);
        cind = job->job;
        val = wctlp_getcoef(lp, &cind, &vind, &cval);
        CCcheck_val_2(val, "Failed wctlp_getcoef");
        cval += 1.0;
        set->nb[job->job] += 1;
        val = wctlp_chgcoef(lp, 1, &cind, &vind, &cval);
        CCcheck_val_2(val, "Failed wctlp_chgcoef");
    }

    cind = njobs;
    cval = 1.0;
    val = wctlp_chgcoef(lp, 1, &cind, &vind, &cval);
    CCcheck_val_2(val, "Failed wctlp_chgcoef");

CLEAN:
    return val;
}

void g_add_col_to_lp(gpointer data, gpointer user_data) {
    scheduleset *tmp = (scheduleset *)data;
    wctdata *    pd = (wctdata *)user_data;
    addColToLP(tmp, pd);
}

int build_lp(wctdata *pd, int construct) {
    int         val = 0;
    wctproblem *problem = pd->problem;
    int         njobs = problem->njobs;
    int         nmachines = problem->nmachines;
    int         nb_row;
    int *       covered = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(covered, "Failed to allocate memory to covered");
    fill_int(covered, njobs, 0);
    val = wctlp_init(&(pd->LP), NULL);
    CCcheck_val_2(val, "wctlp_init failed");

    for (int i = 0; i < njobs; i++) {
        val = wctlp_addrow(pd->LP, 0, (int *)NULL, (double *)NULL, wctlp_GREATER_EQUAL,
                           1.0, (char *)NULL);
        CCcheck_val_2(val, "Failed wctlp_addrow");
    }

    val = wctlp_addrow(pd->LP, 0, (int *)NULL, (double *)NULL, wctlp_LESS_EQUAL,
                       (double)nmachines, (char *)NULL);

    wctlp_get_nb_rows(pd->LP, &nb_row);

    /** add columns from localColPool */
    g_ptr_array_foreach(pd->localColPool, g_add_col_to_lp, pd);

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
    //             CCcheck_NULL_2(pd->newsets[0].C, "Failed to allocate
    //             memory"); pd->newsets[0].table =
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
    //             val = wctlp_addcol(pd->LP, 2, pd->newsets[0].members,
    //             pd->coef,
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
    fill_dbl(pd->pi_sep, nb_row, 0.0);
    pd->subgradient_in = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->subgradient_in, "Failed to allocate memory");
    pd->subgradient = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->subgradient, "Failed to allocate memory");
    pd->x_e = CC_SAFE_MALLOC(2*get_datasize(pd->solver), double);
    CCcheck_NULL_2(pd->x_e, "Failed to allocate memory x_e");
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
    return val;
}

int get_solution_lp_lowerbound(wctdata *pd) {
    int          val = 0;
    int          nbcols;
    scheduleset *tmp;
    Job *        tmp_j;

    val = wctlp_get_nb_cols(pd->LP, &nbcols);
    CCcheck_val_2(val, "Failed in wctlp_get_nb_cols");
    pd->x = CC_SAFE_REALLOC(pd->x, nbcols, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory");
    assert(nbcols == pd->localColPool->len);
    wctlp_x(pd->LP, pd->x, 0);

    for (unsigned i = 0; i < nbcols; ++i) {
        tmp = ((scheduleset *)g_ptr_array_index(pd->localColPool, i));
        if (pd->x[i]) {
            printf("%f: ", pd->x[i]);
            g_ptr_array_foreach(tmp->jobs, g_print_machine, NULL);
            printf("\n");
            for (unsigned j = 0; j < tmp->jobs->len; ++j) {
                tmp_j = (Job *)g_ptr_array_index(tmp->jobs, j);
                printf("%d (%d) ", tmp->nb[tmp_j->job], tmp_j->processingime);
            }
            printf("\n");
        }
    }

CLEAN:
    return val;
}
