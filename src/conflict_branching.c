#include <wct.h>

/** help functions for conflict branching */
static int create_same_conflict(Problem*   problem,
                                NodeData*  parent_pd,
                                NodeData** child,
                                Job*       v1,
                                Job*       v2);
static int create_differ_conflict(Problem*   problem,
                                  NodeData*  parent_pd,
                                  NodeData** child,
                                  Job*       v1,
                                  Job*       v2);

static int transfer_same_cclasses(NodeData*  pd,
                                  GPtrArray* colPool,
                                  Job*       v1,
                                  Job*       v2) {
    int          val = 0;
    ScheduleSet* tmp;
    Job*         tmp_j;
    /* Transfer independent sets: */
    g_ptr_array_set_size(pd->localColPool, colPool->len);

    for (guint i = 0; i < colPool->len; ++i) {
        ScheduleSet* it = (ScheduleSet*)g_ptr_array_index(colPool, i);
        int          construct = 1;
        gboolean     v1_in;
        gboolean     v2_in;
        // v1_in = g_hash_table_contains(it->table, v1);
        // v2_in = g_hash_table_contains(it->table, v2);

        if ((v1_in == 1 && v2_in == 0) || (v1_in == 0 && v2_in == 1)) {
            construct = 0;
        } else {
            tmp = scheduleset_alloc(pd->nb_jobs);
            g_ptr_array_add(pd->localColPool, tmp);
        }

        for (guint j = 0; j < it->job_list->len && construct; ++j) {
            tmp_j = (Job*)g_ptr_array_index(it->job_list, j);
            tmp->total_processing_time += tmp_j->processing_time;
            // g_hash_table_insert(tmp->table, tmp_j, NULL);
            g_ptr_array_add(tmp->job_list, tmp_j);
            tmp->total_weighted_completion_time +=
                value_Fj(tmp->total_processing_time, tmp_j);
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET SAME ");

            g_ptr_array_foreach(it->job_list, g_print_job, NULL);

            printf("\n");
            printf("TRANS SET SAME");

            g_ptr_array_foreach(tmp->job_list, g_print_job, NULL);

            printf("\n");
        }
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
/* END Transfer independent sets: */
CLEAN:
    return val;
}

static int create_same_conflict(Problem*   problem,
                                NodeData*  parent_pd,
                                NodeData** child,
                                Job*       v1,
                                Job*       v2) {
    int val = 0;
    // Parms *parms = &(problem->parms);
    NodeData*   pd = CC_SAFE_MALLOC(1, NodeData);
    Statistics* statistics = parent_pd->stat;
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    nodedata_init(pd, problem);
    /** Init B&B data */
    pd->parent = parent_pd;
    pd->depth = parent_pd->depth + 1;
    // parent_pd->nb_same = 1;
    // parent_pd->same_children = pd;
    pd->v1 = v1;
    pd->v2 = v2;
    /** Init jobs data */
    pd->nb_jobs = parent_pd->nb_jobs;
    pd->nb_machines = parent_pd->nb_machines;
    pd->jobarray = parent_pd->jobarray;
    pd->edge_count_same = parent_pd->edge_count_same + 1;
    pd->edge_count_differ = parent_pd->edge_count_differ;
    /** Init lower bound and upper bound */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->parent = parent_pd;
    // pd->debugcolors = parent_pd->debugcolors;
    // pd->ndebugcolors = parent_pd->ndebugcolors;
    pd->stat = parent_pd->stat;

    /* Construction of solver*/
    if (pd->parent) {
        CCutil_start_resume_time(&(statistics->tot_build_dd));
        // pd->solver = copySolver(pd->parent->solver);
        // add_one_conflict(pd->solver, parms, pd->v1, pd->v2, 1);

        if (pd->nb_jobs != get_num_layers(pd->solver)) {
            pd->status = infeasible;

            if (pd->solver) {
                freeSolver(pd->solver);
                pd->solver = (PricerSolver*)NULL;
            }

            *child = pd;
            CCutil_suspend_timer(&(statistics->tot_build_dd));
            goto CLEAN;
        }

        // init_tables(pd->solver);
        CCutil_suspend_timer(&(statistics->tot_build_dd));
    }

    val = transfer_same_cclasses(pd, parent_pd->localColPool, v1, v2);
    CCcheck_val_2(val, "Failed in transfer_same_cclasses");
    *child = pd;
CLEAN:

    if (val) {
        if (pd) {
            nodedata_free(pd);
            free(pd);
        }
    }

    return val;
}

static int create_differ_conflict(Problem*   problem,
                                  NodeData*  parent_pd,
                                  NodeData** child,
                                  Job*       v1,
                                  Job*       v2) {
    int         val = 0;
    int         i;
    int         nb_cols;
    NodeData*   pd = CC_SAFE_MALLOC(1, NodeData);
    Statistics* statistics = parent_pd->stat;
    // Parms *parms = &(problem->parms);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    ScheduleSet *it, *tmp;
    Job*         tmp_j;
    nodedata_init(pd, problem);
    /** Init B&B data */
    pd->parent = parent_pd;
    pd->depth = parent_pd->depth + 1;
    // parent_pd->nb_diff         += 1;
    // parent_pd->diff_children = pd;
    pd->v1 = v1;
    pd->v2 = v2;
    /** Init jobs data */
    pd->nb_jobs = parent_pd->nb_jobs;
    pd->nb_machines = parent_pd->nb_machines;
    pd->jobarray = parent_pd->jobarray;
    /** Init lower bound and upper bound of node */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    /* Create  graph with extra edge (v1,v2) */
    pd->edge_count_differ = parent_pd->edge_count_differ + 1;
    pd->edge_count_same = parent_pd->edge_count_same;

    /* Construction of solver*/
    if (pd->parent) {
        CCutil_start_resume_time(&(statistics->tot_build_dd));
        // pd->solver = copySolver(pd->parent->solver);
        // add_one_conflict(pd->solver, parms, pd->v1, pd->v2, 0);

        if (pd->nb_jobs != get_num_layers(pd->solver)) {
            pd->status = infeasible;

            if (pd->solver) {
                freeSolver(pd->solver);
                pd->solver = (PricerSolver*)NULL;
            }

            *child = pd;
            CCutil_suspend_timer(&(statistics->tot_build_dd));
            goto CLEAN;
        }

        // init_tables(pd->solver);
        CCutil_suspend_timer(&(statistics->tot_build_dd));
    }

    /* Transfer independent sets by removing v2 if both v1 and v2 are currently
     * contained: */
    nb_cols = parent_pd->localColPool->len;
    pd->localColPool = g_ptr_array_sized_new(nb_cols);

    for (i = 0; i < nb_cols; ++i) {
        it = (ScheduleSet*)g_ptr_array_index(parent_pd->localColPool, i);
        // gboolean v1_in = g_hash_table_contains(it->table, v1);
        // gboolean v2_in = g_hash_table_contains(it->table, v2);
        gboolean v1_in, v2_in;
        int      construct = (v1_in && v2_in) ? 0 : 1;

        if (construct) {
            tmp = scheduleset_alloc(pd->nb_jobs);
            CCcheck_NULL_3(tmp, "Failed to allocate memory");
            g_ptr_array_add(pd->localColPool, tmp);

            for (guint j = 0; j < it->job_list->len; j++) {
                tmp_j = (Job*)g_ptr_array_index(it->job_list, j);
                tmp->total_processing_time += tmp_j->processing_time;
                // g_hash_table_insert(tmp->table, tmp_j, NULL);
                g_ptr_array_add(tmp->job_list, tmp_j);
                tmp->total_weighted_completion_time +=
                    value_Fj(tmp->total_processing_time, tmp_j);
            }
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET DIFFER");

            g_ptr_array_foreach(it->job_list, g_print_machine, NULL);

            printf("\n");
            printf("TRANS SET DIFFER");

            g_ptr_array_foreach(tmp->job_list, g_print_machine, NULL);

            printf("\n");
        }

        CCcheck_val_2(val, "Illegal colorset created in create_differ\n!");
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
    *child = pd;
CLEAN:

    if (val) {
        if (pd) {
            nodedata_free(pd);
            free(pd);
        }
    }

    return val;
}

static int find_strongest_children_conflict(int*           strongest_v1,
                                            int*           strongest_v2,
                                            NodeData*      pd,
                                            Problem*       problem,
                                            HeapContainer* cand_heap,
                                            int*           nodepair_refs,
                                            double*        nodepair_weights) {
    int    val = 0;
    int    max_non_improving_branches = 4; /* pd->nb_jobs / 100 + 1; */
    int    remaining_branches = max_non_improving_branches;
    double strongest_dbl_lb = -115648465146;
    int*   min_nodepair;
    Parms* parms = &(problem->parms);
    *strongest_v1 = -1;
    *strongest_v2 = -1;

    switch (parms->strong_branching) {
        case use_strong_branching:
            while ((min_nodepair = (int*)heapcontainer_min(cand_heap)) &&
                   (remaining_branches--)) {
                int    v1 = -1, v2 = -1;
                double dbl_child_lb;
                inodepair_ref_key(&v1, &v2,
                                  (int)(min_nodepair - nodepair_refs));
                assert(v1 < v2);
                NodeData* same_children = (NodeData*)NULL;
                NodeData* diff_children = (NodeData*)NULL;

                if (dbg_lvl() > 0) {
                    printf(
                        "Creating branches for v1 = %d, v2 = %d (node-pair "
                        "weight %f)\n",
                        v1, v2,
                        nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
                }

                /* Create DIFFER and SAME */
                val = create_same_conflict(
                    problem, pd, &(same_children),
                    (Job*)g_ptr_array_index(problem->g_job_array, v1),
                    (Job*)g_ptr_array_index(problem->g_job_array, v2));
                CCcheck_val_2(val, "Failed in create_same");

                if (same_children->status != infeasible) {
                    same_children->maxiterations = 20;
                    compute_lower_bound(same_children);
                }

                val = create_differ_conflict(
                    problem, pd, &(diff_children),
                    (Job*)g_ptr_array_index(problem->g_job_array, v1),
                    (Job*)g_ptr_array_index(problem->g_job_array, v2));
                CCcheck_val_2(val, "Failed in create_differ");

                if (diff_children->status != infeasible) {
                    diff_children->maxiterations = 20;
                    compute_lower_bound(diff_children);
                }

                dbl_child_lb = (same_children->LP_lower_bound <
                                diff_children->LP_lower_bound)
                                   ? same_children->LP_lower_bound
                                   : diff_children->LP_lower_bound;

                if (dbl_child_lb > strongest_dbl_lb) {
                    strongest_dbl_lb = dbl_child_lb;
                    *strongest_v1 = v1;
                    *strongest_v2 = v2;

                    // remaining_branches = max_non_improving_branches;

                    if (pd->same_children) {
                        nodedata_free(pd->same_children);
                        free(pd->same_children);
                        pd->nb_same = 0;
                    }

                    same_children->maxiterations = 1000000;
                    pd->same_children = same_children;
                    pd->nb_same = 1;

                    if (pd->diff_children) {
                        nodedata_free(pd->diff_children);
                        free(pd->diff_children);
                        pd->nb_diff = 0;
                    }

                    diff_children->maxiterations = 1000000;
                    pd->diff_children = diff_children;
                    pd->nb_diff = 1;
                } else {
                    nodedata_free(same_children);
                    free(same_children);
                    nodedata_free(diff_children);
                    free(diff_children);
                }

                if (dbg_lvl() > 1) {
                    printf(
                        "Found child bound of %f for v1 = %d, v2 = %d, "
                        "nodepair_weight "
                        "= %f .\n",
                        dbl_child_lb, v1, v2,
                        nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
                }
            }

            if (dbg_lvl() > 0) {
                int nodepair_ref =
                    nodepair_ref_key(*strongest_v1, *strongest_v2);
                printf(
                    "Found strongest child bound of %f for v1 = %d, "
                    "v2 = %d, nodepair_weight = %f .\n",
                    strongest_dbl_lb, *strongest_v1, *strongest_v2,
                    nodepair_weights[nodepair_ref]);
            }

            break;

        case no_strong_branching:
            min_nodepair = (int*)heapcontainer_min(cand_heap);
            int v1 = -1, v2 = -1;
            inodepair_ref_key(&v1, &v2, (int)(min_nodepair - nodepair_refs));
            assert(v1 < v2);

            if (dbg_lvl() > 0) {
                printf(
                    "Creating branches for v1 = %d, v2 = %d (node-pair weight "
                    "%f)\n",
                    v1, v2,
                    nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
            }

            val = create_same_conflict(
                problem, pd, &(pd->same_children),
                (Job*)g_ptr_array_index(problem->g_job_array, v1),
                (Job*)g_ptr_array_index(problem->g_job_array, v2));
            CCcheck_val_2(val, "Failed in create_same");
            pd->nb_same = 1;
            val = create_differ_conflict(
                problem, pd, &(pd->diff_children),
                (Job*)g_ptr_array_index(problem->g_job_array, v1),
                (Job*)g_ptr_array_index(problem->g_job_array, v2));
            CCcheck_val_2(val, "Failed in create_differ");
            pd->nb_diff = 1;
            break;
    }

CLEAN:
    return val;
}

int create_branches_conflict(NodeData* pd, Problem* problem) {
    int            val = 0;
    int            status;
    int            i;
    int            nb_cols;
    double*        x = (double*)NULL;
    int            strongest_v1 = -1, strongest_v2 = -1;
    int*           nodepair_refs = (int*)NULL;
    double*        nodepair_weights = (double*)NULL;
    int            nb_pairs = pd->nb_jobs * (pd->nb_jobs + 1) / 2;
    int*           mf_col = (int*)NULL;
    GList*         branchjobs = (GList*)NULL;
    int*           completion_time = (int*)NULL;
    HeapContainer* heap = (HeapContainer*)NULL;
    val = heapcontainer_init(&heap, nb_pairs);
    CCcheck_val_2(val, "Failed pmcheap_init");
    nodepair_refs = CC_SAFE_MALLOC(nb_pairs, int);
    CCcheck_NULL_2(nodepair_refs, "Failed to allocate memory to nodepair_refs");
    nodepair_weights = CC_SAFE_MALLOC(nb_pairs, double);
    CCcheck_NULL_2(nodepair_weights,
                   "Failed to allocate memory to nodepair_weights");

    for (i = 0; i < nb_pairs; i++) {
        nodepair_refs[i] = -1;
        nodepair_weights[i] = .0;
    }

    mf_col = CC_SAFE_MALLOC(pd->nb_jobs, int);
    CCcheck_NULL_2(mf_col, "Failed to allocate memory to mf_col");

    for (i = 0; i < pd->nb_jobs; i++) {
        mf_col[i] = -1;
    }

    if (!pd->RMP) {
        val = build_rmp(pd);
        CCcheck_val_2(val, "Failed at build_lp");
    }

    if (!pd->localColPool->len) {
        compute_lower_bound(pd);
    }

    lp_interface_get_nb_cols(pd->RMP, &nb_cols);
    assert(pd->localColPool->len == nb_cols);
    x = CC_SAFE_MALLOC(nb_cols, double);
    CCcheck_NULL_2(x, "Failed to allocate memory to x");
    val = lp_interface_optimize(pd->RMP, &status);
    CCcheck_val_2(val, "Failed at lp_interface_optimize");

    if (status == GRB_INFEASIBLE) {
        goto CLEAN;
    }

    val = lp_interface_x(pd->RMP, x, 0);
    CCcheck_val_2(val, "Failed at lp_interface_x");
    CC_IFFREE(pd->lambda, double);
    pd->lambda = CC_SAFE_MALLOC(nb_cols, double);
    CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
    memcpy(pd->lambda, x, nb_cols * sizeof(double));
    val = insert_frac_pairs_into_heap(pd, nodepair_refs, nodepair_weights,
                                      nb_pairs, heap);
    CCcheck_val_2(val, "Failed in insert_frac_pairs_into_heap");

    if (heapcontainer_size(heap) == 0) {
        printf("LP returned integral solution\n");
        val = grab_integer_solution(pd, x, lp_int_tolerance());
        CCcheck_val_2(val, "Failed in grab_int_sol");
        assert(pd->status == finished);
        goto CLEAN;
    }

    if (dbg_lvl() > 1) {
        printf("Collected %d branching candidates.\n",
               heapcontainer_size(heap));
    }

    val = find_strongest_children_conflict(&strongest_v1, &strongest_v2, pd,
                                           problem, heap, nodepair_refs,
                                           nodepair_weights);
    CCcheck_val_2(val, "Failed in find_strongest_children");
    // val = create_same_conflict(problem, pd, strongest_v1, strongest_v2);
    // CCcheck_val(val, "Failed in create_same");
    val =
        set_id_and_name(pd->same_children, problem->nb_data_nodes++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");

    if (pd->same_children->status != infeasible) {
        val = compute_lower_bound(pd->same_children);
        CCcheck_val_2(val, "Failed in compute_lower_bound");
    }

    // val = create_differ_conflict(problem, pd, strongest_v1, strongest_v2);
    // CCcheck_val_2(val, "Failed in create_differ");
    val =
        set_id_and_name(pd->diff_children, problem->nb_data_nodes++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");

    if (pd->diff_children->status != infeasible) {
        val = compute_lower_bound(pd->diff_children);
        CCcheck_val_2(val, "Failed in compute_lower_bound");
    }

    if (pd->same_children->status != infeasible &&
        pd->diff_children->status != infeasible) {
        if (pd->same_children->LP_lower_bound <=
            pd->diff_children->LP_lower_bound) {
            pd->same_children->choose = 1;
        } else {
            pd->diff_children->choose = 1;
        }
    }

    free_elist(pd->same_children, &(problem->parms));
    free_elist(pd->diff_children, &(problem->parms));
CLEAN:
    lp_node_data_free(pd);
    free_elist(pd, &(problem->parms));

    if (heap) {
        heapcontainer_free(heap);
        heap = (HeapContainer*)NULL;
    }

    CC_IFFREE(x, double);
    CC_IFFREE(mf_col, int);
    CC_IFFREE(nodepair_refs, int);
    CC_IFFREE(nodepair_weights, double);
    CC_IFFREE(completion_time, int);
    g_list_free(branchjobs);
    return val;
}
