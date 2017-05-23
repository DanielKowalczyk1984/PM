#include <heap.h>
#include <wct.h>

/** help functions for wide branching */
static int create_same_wide(wctproblem *problem,
                            wctdata *   parent_pd,
                            int *       wide_v1,
                            int *       wide_v2,
                            int         nbwide);
static int create_diff_wide(wctproblem *problem,
                            wctdata *   parent_pd,
                            int         v1,
                            int         v2);
static int collect_diff_child_wide(wctdata *cd);
static int collect_same_child_wide(wctdata *cd);
static int remove_finished_subtree_wide(wctdata *child);

// extern inline void inodepair_ref_key(int *v1, int *v2, int index);
// extern inline int nodepair_ref_key(int v1, int v2);

static int collect_same_child_wide(wctdata *cd) {
    int rval = 0;
    int c;

    for (c = 0; c < cd->nsame; ++c) {
        if (cd->same_children_wide[c]->nbbest &&
            (!cd->nbbest ||
             cd->same_children_wide[c]->upper_bound <= cd->upper_bound)) {
            if (cd->nbbest) {
                schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->same_children_wide[c]->upper_bound;
            cd->nbbest = cd->same_children_wide[c]->nbbest;
            cd->same_children_wide[c]->nbbest = 0;
            cd->bestcolors = cd->same_children_wide[c]->bestcolors;
            cd->same_children_wide[c]->bestcolors = (scheduleset *)NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_diff_child_wide(wctdata *cd) {
    int rval = 0;
    int c;

    for (c = 0; c < cd->ndiff; ++c) {
        if (cd->diff_children_wide[c]->nbbest &&
            (!cd->nbbest ||
             cd->diff_children_wide[c]->upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct =
                cd->diff_children_wide[c]->besttotwct;
            cd->nbbest = cd->diff_children_wide[c]->nbbest;
            cd->diff_children_wide[c]->nbbest = 0;
            cd->bestcolors = cd->diff_children_wide[c]->bestcolors;
            cd->diff_children_wide[c]->bestcolors = (scheduleset *)NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int remove_finished_subtree_wide(wctdata *child) {
    int      val = 0;
    int      i;
    wctdata *cd = (wctdata *)child;
    int      all_same_finished = 1;
    int      all_diff_finished = 1;

    while (cd) {
        for (i = 0; i < cd->nsame; ++i) {
            if (cd->same_children_wide[i]->status < finished) {
                all_same_finished = 0;
                break;
            }
        }

        if (cd->nsame && all_same_finished) {
            val = collect_same_child_wide(cd);
            CCcheck_val_2(val, "Failed in collect_same_children");

            for (i = 0; i < cd->nsame; ++i) {
                wctdata_free(cd->same_children_wide[i]);
                CC_IFFREE(cd->same_children_wide[i], wctdata);
            }

            free(cd->same_children_wide);
            cd->same_children_wide = (wctdata **)NULL;
            cd->nsame = 0;
        }

        for (i = 0; i < cd->ndiff; ++i) {
            if (cd->diff_children_wide[i]->status < finished) {
                all_diff_finished = 0;
                break;
            }
        }

        if (cd->ndiff && all_diff_finished) {
            val = collect_diff_child_wide(cd);
            CCcheck_val_2(val, "Failed in collect_diff_children");

            for (i = 0; i < cd->ndiff; ++i) {
                wctdata_free(cd->diff_children_wide[i]);
                CC_IFFREE(cd->diff_children_wide[i], wctdata);
            }

            CC_IFFREE(cd->diff_children_wide, wctdata *);
            cd->ndiff = 0;
        }

        if (!cd->same_children && !cd->diff_children) {
            cd->status = finished;
            CCcheck_val_2(val, "Failed to write_wctdata");
            cd = cd->parent;
        } else {
            cd = (wctdata *)NULL;
        }
    }

CLEAN:
    return val;
}

static int transfer_same_cclasses_wide(wctdata *          pd,
                                       const scheduleset *parent_cclasses,
                                       int                parent_ccount,
                                       int *              v1_wide,
                                       int *              v2_wide) {
    int val = 0;
    int i;
    /* Transfer independent sets: */
    pd->gallocated = pd->ccount = parent_ccount;
    pd->cclasses = CC_SAFE_MALLOC(pd->gallocated, scheduleset);
    CCcheck_NULL_2(pd->cclasses, "Failed to allocate memory");
    pd->ccount = 0;

    for (i = 0; i < parent_ccount; ++i) {
        int j;
        int construct = 1;

        for (j = 0; j < pd->nb_wide && construct; j++) {
            gboolean v1_in = g_hash_table_contains(parent_cclasses[i].table,
                                                   GINT_TO_POINTER(v1_wide[i]));
            gboolean v2_in = g_hash_table_contains(parent_cclasses[i].table,
                                                   GINT_TO_POINTER(v2_wide[i]));

            if ((v1_in == 1 && v2_in == 0) || (v1_in == 0 && v2_in == 1)) {
                construct = 0;
                break;
            }
        }

        if (construct) {
            scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members =
                CC_SAFE_MALLOC(parent_cclasses[i].count + 1, int);
            pd->cclasses[pd->ccount].C =
                CC_SAFE_MALLOC(parent_cclasses[i].count, int);
            pd->cclasses[pd->ccount].table =
                g_hash_table_new(g_direct_hash, g_direct_equal);
            pd->cclasses[pd->ccount].count = 0;
            pd->cclasses[pd->ccount].age = 0;
            pd->cclasses[pd->ccount].totweight = 0;
            pd->cclasses[pd->ccount].totwct = 0;
        } else {
            continue;
        }

        for (j = 0; j < parent_cclasses[i].count; ++j) {
            pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)] =
                parent_cclasses[i].members[j];
            pd->cclasses[pd->ccount].totweight +=
                pd->jobarray[parent_cclasses[i].members[j]].processingime;
            pd->cclasses[pd->ccount].C[(pd->cclasses[pd->ccount].count)] =
                pd->cclasses[pd->ccount].totweight;
            g_hash_table_insert(
                pd->cclasses[pd->ccount].table,
                GINT_TO_POINTER(pd->cclasses[pd->ccount]
                                    .members[(pd->cclasses[pd->ccount].count)]),
                pd->cclasses[pd->ccount].C + (pd->cclasses[pd->ccount].count));
            pd->cclasses[pd->ccount].totwct +=
                pd->jobarray[parent_cclasses[i].members[j]].weight *
                pd->cclasses[pd->ccount].totweight;
            pd->cclasses[pd->ccount].count++;
            /* else 'parent_cclasses[i].members[j] == v2' and we skip it*/
        }

        pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] =
            pd->njobs;
        pd->ccount++;

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET SAME ");

            for (j = 0; j < parent_cclasses[i].count; ++j) {
                printf(" %d", parent_cclasses[i].members[j]);
            }

            printf("\n");
            printf("TRANS SET SAME");

            for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
            }

            printf("\n");
        }
    }

    for (i = pd->ccount; i < pd->gallocated; i++) {
        scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
/* END Transfer independent sets: */
CLEAN:
    return val;
}

static int create_diff_wide(wctproblem *problem,
                            wctdata *   parent_pd,
                            int         v1,
                            int         v2) {
    int       val = 0;
    int       i;
    wctdata * pd = CC_SAFE_MALLOC(1, wctdata);
    wctparms *parms = &(problem->parms);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd, problem);
    /** Init B&B data */
    pd->parent = parent_pd;
    pd->depth = parent_pd->depth + 1;
    parent_pd->diff_children_wide[parent_pd->ndiff++] = pd;
    pd->v1 = v1;
    pd->v2 = v2;
    /** Init jobs data */
    pd->njobs = parent_pd->njobs;
    pd->nmachines = parent_pd->nmachines;
    pd->jobarray = parent_pd->jobarray;
    /** Init lower bound and upper bound of node */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_pd->dbl_safe_lower_bound;
    /* Create  graph with extra edge (v1,v2) */
    pd->ecount_differ = parent_pd->ecount_differ + 1;
    pd->ecount_same = parent_pd->ecount_same;
    // pd->elist_differ  = CC_SAFE_MALLOC(2 * pd->ecount_differ, int);
    // CCcheck_NULL_2(pd->elist_differ, "Failed to allocate pd->elist");
    // if (parent_pd->ecount_differ > 0) {
    //     memcpy(pd->elist_differ, parent_pd->elist_differ, 2 *
    //     parent_pd->ecount_differ * sizeof(int));
    // }
    // pd->elist_differ[ 2 * (pd->ecount_differ - 1)] = v1;
    // pd->elist_differ[ 2 * (pd->ecount_differ - 1) + 1] = v2;
    // pd->debugcolors = parent_pd->debugcolors;
    // pd->ndebugcolors = parent_pd->ndebugcolors;
    // /** Copy same list */
    // if (parent_pd->ecount_same > 0) {
    //     pd->elist_same = CC_SAFE_MALLOC(2 * pd->ecount_same, int);
    //     memcpy(pd->elist_same, parent_pd->elist_same, 2 * pd->ecount_same *
    //     sizeof(int));
    // }
    // /* END: Create  graph with extra edge (v1,v2) */
    // if (dbg_lvl() > 1) {
    //     printf("create_differ created following graph:\n");
    //     adjGraph_print(pd->ecount_differ, pd->elist_differ);
    // }

    /* Construction of solver*/
    if (pd->parent &&
        (parms->solver == bdd_solver || parms->solver == zdd_solver)) {
        CCutil_start_resume_time(&(problem->tot_build_dd));
        pd->solver = copySolver(pd->parent->solver);
        add_one_conflict(pd->solver, parms, pd->v1, pd->v2, 0);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case DP_solver:
                break;
        }

        init_tables(pd->solver);
        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    /* Transfer independent sets by removing v2 if both v1 and v2 are currently
     * contained: */
    pd->gallocated = parent_pd->ccount;
    pd->cclasses = CC_SAFE_MALLOC(pd->gallocated, scheduleset);
    pd->ccount = 0;

    for (i = 0; i < parent_pd->ccount; ++i) {
        int j;
        int v1_found = 0;
        int construct = 1;

        for (j = 0; j < parent_pd->cclasses[i].count; ++j) {
            int current_elm = parent_pd->cclasses[i].members[j];

            if (current_elm == v1) {
                v1_found = 1;
            }

            if (current_elm == v2) {
                if (v1_found) {
                    construct = 0;
                }
            }
        }

        if (construct) {
            scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members =
                CC_SAFE_MALLOC(parent_pd->cclasses[i].count + 1, int);
            CCcheck_NULL_2(pd->cclasses[pd->ccount].members,
                           "Failed to allocate pd->cclasses[i].members");
            pd->cclasses[pd->ccount].C =
                CC_SAFE_MALLOC(parent_pd->cclasses[i].count, int);
            CCcheck_NULL_2(pd->cclasses[pd->ccount].members,
                           "Failed to allocate memory");
            pd->cclasses[pd->ccount].table =
                g_hash_table_new(g_direct_hash, g_direct_equal);

            for (j = 0; j < parent_pd->cclasses[i].count; j++) {
                pd->cclasses[pd->ccount]
                    .members[pd->cclasses[pd->ccount].count] =
                    parent_pd->cclasses[i].members[j];
                pd->cclasses[pd->ccount].totweight +=
                    parent_pd->jobarray[parent_pd->cclasses[i].members[j]]
                        .processingime;
                pd->cclasses[pd->ccount].C[pd->cclasses[pd->ccount].count] =
                    pd->cclasses[pd->ccount].totweight;
                g_hash_table_insert(
                    pd->cclasses[pd->ccount].table,
                    GINT_TO_POINTER(
                        pd->cclasses[pd->ccount]
                            .members[pd->cclasses[pd->ccount].count]),
                    pd->cclasses[pd->ccount].C +
                        pd->cclasses[pd->ccount].count);
                pd->cclasses[pd->ccount].totwct +=
                    parent_pd->jobarray[parent_pd->cclasses[i].members[j]]
                        .weight *
                    pd->cclasses[pd->ccount].totweight;
                (pd->cclasses[pd->ccount].count)++;
            }

            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] =
                pd->njobs;
            pd->ccount++;
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET DIFFER");

            for (j = 0; j < parent_pd->cclasses[i].count; ++j) {
                printf(" %d", parent_pd->cclasses[i].members[j]);
            }

            printf("\n");
            printf("TRANS SET DIFFER");

            for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
            }

            printf("\n");
        }

        CCcheck_val_2(val, "Illegal colorset created in create_differ\n!");
    }

    for (i = pd->ccount; i < pd->gallocated; i++) {
        scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_pd->diff_children = (wctdata *)NULL;
    }

    return val;
}

static int create_same_wide(wctproblem *problem,
                            wctdata *   parent_pd,
                            int *       v1_wide,
                            int *       v2_wide,
                            int         nb_wide) {
    int       val = 0;
    wctparms *parms = &(problem->parms);
    wctdata * pd = CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd, problem);
    /** Init B&B data */
    pd->depth = parent_pd->depth + 1;
    parent_pd->same_children_wide[parent_pd->nsame++] = pd;
    /** Init jobs data */
    pd->njobs = parent_pd->njobs;
    pd->nmachines = parent_pd->nmachines;
    pd->jobarray = parent_pd->jobarray;
    pd->ecount_same = parent_pd->ecount_same + 1;
    pd->ecount_differ = parent_pd->ecount_differ;
    /** Init lower bound and upper bound */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_pd->dbl_safe_lower_bound;
    pd->parent = parent_pd;
    pd->debugcolors = parent_pd->debugcolors;
    pd->ndebugcolors = parent_pd->ndebugcolors;
    pd->nb_wide = nb_wide;
    pd->v1_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(pd->v1_wide, "Failed to allocate memory");
    pd->v2_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(pd->v2_wide, "Failed to allocate memory");
    memcpy(pd->v1_wide, v1_wide, sizeof(int) * nb_wide);
    memcpy(pd->v2_wide, v2_wide, sizeof(int) * nb_wide);
    pd->elist_same = CC_SAFE_MALLOC(2 * pd->nb_wide, int);
    CCcheck_NULL_2(pd->elist_same, "Failed to allocate memory");

    for (int i = 0; i < nb_wide; i++) {
        pd->elist_same[2 * i] = v1_wide[i];
        pd->elist_same[2 * i + 1] = v2_wide[i];
    }

    pd->ecount_same += nb_wide;

    /* Construction of solver*/
    if (pd->parent &&
        (parms->solver == bdd_solver || parms->solver == zdd_solver)) {
        CCutil_start_resume_time(&(problem->tot_build_dd));
        pd->solver = copySolver(pd->parent->solver);
        add_conflict_constraints(pd->solver, &(problem->parms), pd->elist_same,
                                 pd->nb_wide, NULL, 0);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case DP_solver:
                break;
        }

        init_tables(pd->solver);
        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    val = transfer_same_cclasses_wide(
        pd, parent_pd->cclasses, parent_pd->ccount, pd->v1_wide, pd->v2_wide);
    CCcheck_val_2(val, "Failed in transfer_same_cclasses");
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_pd->same_children = (wctdata *)NULL;
    }

    CC_IFFREE(pd->elist_same, int);
    pd->ecount_same = 0;
    return val;
}

int create_branches_wide(wctdata *pd, wctproblem *problem) {
    int       val = 0;
    int       status;
    int       i;
    int *     v1_wide = (int *)NULL;
    int *     v2_wide = (int *)NULL;
    int *     min_nodepair;
    int       nb_wide;
    wctparms *parms = &(problem->parms);
    int *     nodepair_refs = (int *)NULL;
    double *  nodepair_weights = (double *)NULL;
    int       npairs = pd->njobs * (pd->njobs + 1) / 2;
    int *     mf_col = (int *)NULL;
    pmcheap * heap = (pmcheap *)NULL;
    val = pmcheap_init(&heap, npairs);
    CCcheck_val_2(val, "Failed pmcheap_init");
    nodepair_refs = CC_SAFE_MALLOC(npairs, int);
    CCcheck_NULL_2(nodepair_refs, "Failed to allocate memory to nodepair_refs");
    nodepair_weights = CC_SAFE_MALLOC(npairs, double);
    CCcheck_NULL_2(nodepair_weights,
                   "Failed to allocate memory to nodepair_weights");

    for (i = 0; i < npairs; i++) {
        nodepair_refs[i] = -1;
        nodepair_weights[i] = .0;
    }

    if (!pd->LP) {
        val = build_lp(pd, parms->construct);
        CCcheck_val_2(val, "Failed at build_lp");
    }

    if (!pd->ccount) {
        compute_lower_bound(problem, pd);
    }

    assert(pd->ccount != 0);
    val = wctlp_optimize(pd->LP, &status);
    CCcheck_val_2(val, "Failed at wctlp_optimize");

    if (status == GRB_INFEASIBLE) {
        goto CLEAN;
    }

    CC_IFFREE(pd->x, double);
    pd->x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    val = wctlp_x(pd->LP, pd->x, 0);
    CCcheck_val_2(val, "Failed at wctlp_x");
    val = insert_frac_pairs_into_heap(pd, pd->x, nodepair_refs,
                                      nodepair_weights, npairs, heap);
    CCcheck_val_2(val, "Failed in insert_frac_pairs_into_heap");

    if (pmcheap_size(heap) == 0) {
        printf("LP returned integral solution\n");
        val = grab_int_sol(pd, pd->x, lp_int_tolerance());
        CCcheck_val_2(val, "Failed in grab_int_sol");
        assert(pd->status == finished);
        goto CLEAN;
    }

    if (dbg_lvl() > 1) {
        printf("Collected %d branching candidates.\n", pmcheap_size(heap));
    }

    nb_wide = CC_MIN(pmcheap_size(heap), 5);
    pd->same_children_wide = CC_SAFE_MALLOC(1, wctdata *);
    CCcheck_NULL_2(pd->same_children_wide, "Failed to allocate memory");
    pd->diff_children_wide = CC_SAFE_MALLOC(nb_wide, wctdata *);
    CCcheck_NULL_2(pd->diff_children_wide, "Failed to allocate memory");
    v1_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(v1_wide, "Failed to allocate");
    v2_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(v2_wide, "Failed to allocate");

    for (i = 0; i < nb_wide; i++) {
        min_nodepair = (int *)pmcheap_min(heap);
        int v1 = -1, v2 = -1;
        inodepair_ref_key(&v1, &v2, (int)(min_nodepair - nodepair_refs));
        assert(v1 < v2);
        v1_wide[i] = v1;
        v2_wide[i] = v2;

        if (dbg_lvl() > 1) {
            printf(
                "Inserted  v1 = %d and v2 = %d with nodepair_weight = %f .\n",
                v1, v2, nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
        }
    }

    val = create_same_wide(problem, pd, v1_wide, v2_wide, nb_wide);
    CCcheck_val(val, "Failed in create_same");
    /** Set name and id for same children */
    val = compute_lower_bound(problem, *(pd->same_children_wide));
    CCcheck_val_2(val, "Failed in compute_lower_bound");

    for (i = 0; i < nb_wide; ++i) {
        create_diff_wide(problem, pd, v1_wide[i], v2_wide[i]);
        CCcheck_val_2(val, "Failed in create_differ");
        val = set_id_and_name(pd->diff_children_wide[i], problem->nwctdata++,
                              pd->pname);
        CCcheck_val_2(val, "Failed in set_id_and_name");
        val = compute_lower_bound(problem, pd->diff_children_wide[i]);
        CCcheck_val_2(val, "Failed in compute_lower_bound");
        free_elist(pd->diff_children_wide[i], &(problem->parms));
    }

    free_elist(*(pd->same_children_wide), &(problem->parms));
CLEAN:
    lpwctdata_free(pd);
    free_elist(pd, &(problem->parms));

    if (heap) {
        pmcheap_free(heap);
        heap = (pmcheap *)NULL;
    }

    CC_IFFREE(v1_wide, int);
    CC_IFFREE(v2_wide, int);
    CC_IFFREE(mf_col, int);
    CC_IFFREE(nodepair_refs, int);
    CC_IFFREE(nodepair_weights, double);
    return val;
}

int branching_msg_wide(wctdata *pd, wctproblem *problem) {
    BinomialHeap *heap = problem->br_heap_a;

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf(
            "Branching with lb %d (LP %f) at depth %d (id = %d, "
            "time = %f, unprocessed nodes = %u, nbjobs= %d, upper bound = %d, "
            "lower bound = %d, v1 = %d, v2 = %d, nbdiff = %d, nbsame = %d ).\n",
            pd->lower_bound, pd->LP_lower_bound, pd->depth, pd->id,
            problem->tot_cputime.cum_zeit, binomial_heap_num_entries(heap),
            pd->njobs, problem->global_upper_bound, problem->global_lower_bound,
            pd->v1, pd->v2, pd->ecount_differ, pd->ecount_same);
        CCutil_resume_timer(&problem->tot_cputime);
        problem->nb_explored_nodes++;
    }

    return 0;
}

int sequential_branching_wide(wctproblem *problem) {
    int           val = 0;
    wctdata *     pd;
    BinomialHeap *br_heap = problem->br_heap_a;
    wctparms *    parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL WIDE BRANCHING:\n");
    CCutil_suspend_timer(&problem->tot_branch_and_bound);

    while ((pd = (wctdata *)binomial_heap_pop(br_heap)) &&
           problem->tot_branch_and_bound.cum_zeit <
               parms->branching_cpu_limit) {
        CCutil_resume_timer(&problem->tot_branch_and_bound);
        int i;
        pd->upper_bound = problem->global_upper_bound;

        if (pd->lower_bound >= pd->upper_bound ||
            pd->eta_in > pd->upper_bound - 1) {
            skip_wctdata(pd, problem);
            remove_finished_subtree_wide(pd);
        } else {
            branching_msg_wide(pd, problem);
            /** Construct PricerSolver */
            /*val = recover_elist(pd);
            CCcheck_val_2(val, "Failed in recover_elist");*/

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches_wide(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nsame; i++) {
                val = insert_into_branching_heap((pd->same_children_wide[i]),
                                                 problem);
                CCcheck_val_2(val, "Failed in insert_into_branching_heap");
            }

            for (i = 0; i < pd->ndiff; i++) {
                val = insert_into_branching_heap(pd->diff_children_wide[i],
                                                 problem);
                CCcheck_val_2(val, "Faield at insert_into_branching_heap");
            }

            assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                remove_finished_subtree_wide(pd);
            }
        }

        CCutil_suspend_timer(&problem->tot_branch_and_bound);
    }

    CCutil_resume_timer(&problem->tot_branch_and_bound);

    if (pd) {
        printf("Branching timeout of %f second reached\n",
               parms->branching_cpu_limit);
    }

    children_data_free(&problem->root_pd);
CLEAN:
    return val;
}
