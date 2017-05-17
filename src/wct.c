#include <wct.h>

int debug = 0;

/*Information about debug*/
int  dbg_lvl() { return debug; }
void set_dbg_lvl(int dbglvl) { debug = dbglvl; }

/*Functions for initialization of the problem and freeing the problem*/
void wctproblem_init(wctproblem *problem) {
    /** Job data */
    problem->g_job_array = g_ptr_array_new_with_free_func(free);
    problem->opt_sol = (solution *)NULL;
    /** Job summary */
    problem->njobs = 0;
    problem->psum = 0;
    problem->pmax = 0;
    problem->pmin = INT_MAX;
    problem->dmax = INT_MIN;
    problem->dmin = INT_MAX;
    problem->off = 0;
    problem->H_min = 0;
    problem->H_max = INT_MAX;
    problem->e = g_ptr_array_new_with_free_func(g_interval_free);
    /*B&B info*/
    problem->nwctdata = 0;
    problem->global_upper_bound = INT_MAX;
    problem->first_lower_bound = 0;
    problem->global_lower_bound = 0;
    problem->first_upper_bound = INT_MAX;
    problem->status = no_sol;
    problem->rel_error = 1.0;
    problem->nbestschedule = 0;
    problem->bestschedule = (scheduleset *)NULL;
    problem->maxdepth = 0;
    problem->nbinitsets = 0;
    problem->gallocated = 0;
    problem->initsets = (scheduleset *)NULL;
    problem->br_heap_a = (BinomialHeap *)NULL;
    /*data of the problem*/
    wctdata_init(&(problem->root_pd));
    /*parms of the problem*/
    wctparms_init(&(problem->parms));
    /*heap initialization*/
    problem->unexplored_states = g_ptr_array_new();
    problem->non_empty_level_pqs = g_queue_new();
    problem->last_explored = -1;
    problem->found = 0;
    problem->nb_explored_nodes = 0;
    problem->nb_generated_col = 0;
    /*CPU timer initialisation*/
    CCutil_init_timer(&(problem->tot_cputime), "tot_cputime");
    CCutil_init_timer(&(problem->tot_scatter_search), "tot_scatter_search");
    CCutil_init_timer(&(problem->tot_branch_and_bound), "tot_branch_and_bound");
    CCutil_init_timer(&(problem->tot_lb_lp_root), "tot_lb_lp_root");
    CCutil_init_timer(&(problem->tot_build_dd), "tot_build_dd");
    CCutil_init_timer(&(problem->tot_lb_lp), "tot_lb_lp");
    CCutil_init_timer(&(problem->tot_lb), "tot_lb");
    CCutil_init_timer(&(problem->tot_pricing), "tot_pricing");
}

void wctproblem_free(wctproblem *problem) {
    /*free the parameters*/
    wctparms_free(&(problem->parms));
    wctdata_free(&(problem->root_pd));

    /*free the heap*/
    if (problem->br_heap_a != (BinomialHeap *)NULL) {
        binomial_heap_free(problem->br_heap_a);
    }

    for (unsigned int i = 0; i < problem->unexplored_states->len; ++i) {
        if ((problem->unexplored_states->pdata[i]) != (BinomialHeap *)NULL) {
            binomial_heap_free(
                (BinomialHeap *)problem->unexplored_states->pdata[i]);
        }
    }

    g_ptr_array_free(problem->g_job_array, TRUE);
    g_ptr_array_free(problem->e, TRUE);
    g_ptr_array_free(problem->unexplored_states, TRUE);
    g_queue_free(problem->non_empty_level_pqs);
    schedulesets_free(&(problem->initsets), &(problem->gallocated));
    schedulesets_free(&(problem->bestschedule), &(problem->nbestschedule));
    solution_free(&(problem->opt_sol));
}

/*Functions for initialization and free the data*/
void wctdata_init(wctdata *pd) {
    /*Initialization B&B data*/
    pd->id = -1;
    pd->depth = 0;
    pd->status = initialized;
    sprintf(pd->pname, "temporary");
    /*Initialization graph data*/
    pd->njobs = 0;
    pd->orig_node_ids = (int *)NULL;
    // pd->duration = (int *)NULL;
    // pd->weights = (int *)NULL;
    // pd->duetime = (int *) NULL;
    // pd->releasetime = (int *) NULL;
    //pd->jobarray = (Job *)NULL;
    pd->H_max = 0;
    pd->H_min = 0;
    pd->upper_bound = INT_MAX;
    pd->lower_bound = 0;
    pd->dbl_safe_lower_bound = 0.0;
    pd->dbl_est_lower_bound = 0.0;
    pd->lower_scaled_bound = 1;
    pd->kpc_pi_scalef = 1;
    pd->LP_lower_bound = 0.0;
    pd->partial_sol = 0.0;
    pd->rhs = (double *)NULL;
    /*Initialization  of the LP*/
    pd->LP = (wctlp *)NULL;
    pd->x = (double *)NULL;
    pd->coef = (double *)NULL;
    pd->pi = (double *)NULL;
    pd->kpc_pi = (int *)NULL;
    /**init stab data */
    pd->pi_in = (double *)NULL;
    pd->pi_out = (double *)NULL;
    pd->pi_sep = (double *)NULL;
    pd->subgradient = (double *)NULL;
    pd->subgradient_in = (double *)NULL;
    pd->alpha = 0.8;
    pd->update = 1;
    /*Initialization pricing_problem*/
    pd->solver = (PricerSolver *)NULL;
    pd->nnonimprovements = 0;
    /*Initialization of scheduleset*/
    pd->ccount = 0;
    pd->cclasses = (scheduleset *)NULL;
    pd->dzcount = 0;
    pd->gallocated = 0;
    pd->newsets = (scheduleset *)NULL;
    pd->nnewsets = 0;
    pd->bestcolors = (scheduleset *)NULL;
    pd->nbbest = 0;
    pd->debugcolors = (scheduleset *)NULL;
    pd->ndebugcolors = 0;
    pd->opt_track = 0;
    /*Initialization max and retirement age*/
    pd->maxiterations = 1000000;
    pd->retirementage = 1000000;
    /*initialization of branches*/
    pd->branch_job = -1;
    pd->parent = (wctdata *)NULL;
    pd->choose = 0;
    /** ahv branching */
    pd->duetime_child = (wctdata *)NULL;
    pd->nduetime = 0;
    pd->releasetime_child = (wctdata *)NULL;
    pd->nreleasetime = 0;
    pd->branch_job = -1;
    pd->completiontime = 0;
    /** conflict branching */
    pd->elist_same = (int *)NULL;
    pd->ecount_same = 0;
    pd->elist_differ = (int *)NULL;
    pd->ecount_differ = 0;
    pd->same_children = (wctdata *)NULL;
    pd->nsame = 0;
    pd->diff_children = (wctdata *)NULL;
    pd->ndiff = 0;
    pd->v1 = -1;
    pd->v2 = -1;
    /** Wide branching */
    pd->v1_wide = (int *)NULL;
    pd->v2_wide = (int *)NULL;
    pd->nb_wide = 0;
    pd->same_children_wide = (wctdata **)NULL;
    pd->diff_children_wide = (wctdata **)NULL;
}

void lpwctdata_free(wctdata *pd) {
    if (pd->LP) {
        wctlp_free(&(pd->LP));
    }

    CC_IFFREE(pd->coef, double);
    CC_IFFREE(pd->pi, double);
    CC_IFFREE(pd->x, double);
    CC_IFFREE(pd->kpc_pi, int);
    CC_IFFREE(pd->pi_out, double);
    CC_IFFREE(pd->pi_in, double);
    CC_IFFREE(pd->pi_sep, double);
    CC_IFFREE(pd->subgradient, double);
    CC_IFFREE(pd->subgradient_in, double);
    CC_IFFREE(pd->rhs, double);

    if (pd->solver) {
        freeSolver(pd->solver);
        pd->solver = (PricerSolver *)NULL;
    }

    schedulesets_free(&(pd->newsets), &(pd->nnewsets));
    schedulesets_free(&(pd->cclasses), &(pd->gallocated));
    pd->ccount = 0;
}

void children_data_free(wctdata *pd) {
    int i;

    for (i = 0; i < pd->nsame; ++i) {
        wctdata_free(&(pd->same_children[i]));
    }

    for (i = 0; i < pd->nduetime; ++i) {
        wctdata_free(&(pd->duetime_child[i]));
    }

    CC_IFFREE(pd->same_children, wctdata);
    CC_IFFREE(pd->duetime_child, wctdata);

    for (i = 0; i < pd->ndiff; ++i) {
        wctdata_free(&(pd->diff_children[i]));
    }

    for (i = 0; i < pd->nreleasetime; ++i) {
        wctdata_free(&(pd->releasetime_child[i]));
    }

    CC_IFFREE(pd->releasetime_child, wctdata);
    CC_IFFREE(pd->diff_children, wctdata);
    pd->nsame = pd->ndiff = 0;
    pd->nreleasetime = pd->nduetime = 0;
}

void temporary_data_free(wctdata *pd) {
    children_data_free(pd);
    lpwctdata_free(pd);
}

void wctdata_free(wctdata *pd) {
    schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    temporary_data_free(pd);
    CC_IFFREE(pd->elist_same, int);
    CC_IFFREE(pd->elist_differ, int);
    CC_IFFREE(pd->v1_wide, int);
    CC_IFFREE(pd->v2_wide, int);
    CC_IFFREE(pd->orig_node_ids, int);
    CC_IFFREE(pd->v1_wide, int);
    CC_IFFREE(pd->v2_wide, int);
}

int partlist_to_scheduleset(
    partlist *part, int nbpart, int njobs, scheduleset **classes, int *ccount) {
    int val = 0;
    int i, j = 0, k = 0;
    schedulesets_free(classes, ccount);
    *ccount = j;
    scheduleset *temp_sets = CC_SAFE_MALLOC(*ccount, scheduleset);
    CCcheck_NULL_2(temp_sets, "Failed to allocate memory to temp_sets");

    for (i = 0; i < nbpart; i++) {
        if (0) {
            scheduleset_init(temp_sets + k);
            temp_sets[k].totweight = part[i].c;
            temp_sets[k].members = CC_SAFE_MALLOC(temp_sets[k].count + 1, int);
            CCcheck_NULL(temp_sets[k].members,
                         "Failed to allocate memory to members");
            temp_sets[k].totwct = 0;
            temp_sets[k].C = CC_SAFE_MALLOC(temp_sets[k].count, int);
            CCcheck_NULL(temp_sets[k].C, "Failed to allocate memory to C");
            temp_sets[k].table =
                g_hash_table_new(g_direct_hash, g_direct_equal);
            CCcheck_NULL(temp_sets[k].table, "Failed to construct table");
            GList *v = (GList *)NULL;
            j = 0;
            int completion = 0;

            while (v != NULL) {
                Job *job = (Job *)v->data;
                completion += job->processingime;
                temp_sets[k].C[j] = completion;
                temp_sets[k].members[j] = job->job;
                g_hash_table_insert(temp_sets[k].table,
                                    GINT_TO_POINTER(job->job),
                                    temp_sets[k].C + j);
                temp_sets[k].totwct += job->weight * completion;
                v = g_list_next(v);
                j++;
            }

            temp_sets[k].members[j] = njobs;
            // CCutil_int_array_quicksort( temp_sets[k].members,
            // temp_sets[k].count );
            k++;
        }
    }

    *classes = temp_sets;
CLEAN:

    if (val) {
        schedulesets_free(&temp_sets, ccount);
    }

    return val;
}

int set_id_and_name(wctdata *pd, int id, const char *fname) {
    int val = 0;
    int sval = 0;
    CCcheck_NULL_2(pd, "np memory was allocated to pd");
    pd->id = id;
    sval = snprintf(pd->pname, MAX_PNAME_LEN, "%s", fname);

    if (sval < 0 || MAX_PNAME_LEN <= sval) {
        val = 1;
        CCcheck_val(val, "Failed to write pname")
    }

CLEAN:
    return val;
}

static int prefill_heap(wctdata *pd, wctproblem *problem) {
    int val = 0;
    int insert_into_heap = 0;

    if (problem->nwctdata <= pd->id) {
        problem->nwctdata = pd->id + 1;
    }

    if (pd->status < LP_bound_computed) {
        printf("Found a node with LP not computed!\n");
        val = compute_lower_bound(problem, pd);
        CCcheck_val_2(val, "Failed at compute_lower_bound");
        insert_into_heap = 1;
    }

    if (pd->status < finished) {
        int i;

        if (!pd->nsame || !pd->ndiff) {
            insert_into_heap = 1;
        }

        for (i = 0; (!insert_into_heap) && i < pd->nsame; ++i) {
            if (pd->duetime_child[i].status < LP_bound_computed) {
                insert_into_heap = 1;
            }
        }

        for (i = 0; (!insert_into_heap) && i < pd->ndiff; ++i) {
            if (pd->releasetime_child[i].status < LP_bound_computed) {
                insert_into_heap = 1;
            }
        }
    }

    if (insert_into_heap) {
        val = insert_into_branching_heap(pd, problem);
        CCcheck_val_2(val, "Failed in insert_into_branching_heap");
        children_data_free(pd);
    } else {
        int i;

        for (i = 0; i < pd->nsame; ++i) {
            prefill_heap(pd->duetime_child + i, problem);
        }

        for (i = 0; i < pd->ndiff; ++i) {
            prefill_heap(pd->releasetime_child + i, problem);
        }
    }

CLEAN:
    return val;
}

int compute_schedule(wctproblem *problem) {
    int       val = 0;
    wctdata * root_pd = &(problem->root_pd);
    wctparms *parms = &(problem->parms);
    problem->mult_key = 1.0;
    problem->first_upper_bound = problem->global_upper_bound;
    problem->first_lower_bound = problem->global_lower_bound;
    problem->first_rel_error =
        (double)(problem->global_upper_bound - problem->global_lower_bound) /
        ((double)problem->global_lower_bound);
    prune_duplicated_sets(root_pd);
    init_BB_tree(problem);
    print_size_to_csv(problem, root_pd);

    if (root_pd->status >= LP_bound_computed) {
        val = prefill_heap(root_pd, problem);
        CCcheck_val(val, "Failed in prefill_heap");
    } else {
        CCutil_start_timer(&(problem->tot_lb_lp_root));
        val = compute_lower_bound(problem, root_pd);
        CCcheck_val_2(val, "Failed in compute_lower_bound");

        if (root_pd->lower_bound > problem->global_lower_bound) {
            problem->global_lower_bound = root_pd->lower_bound;
            problem->first_lower_bound = root_pd->lower_bound;
            problem->first_rel_error = (double)(problem->first_upper_bound -
                                                problem->first_lower_bound) /
                                       ((double)problem->first_lower_bound);
        }

        problem->parms.construct = 1;
        CCcheck_val_2(val, "Failed in compute_lower_bound");
        problem->nb_generated_col_root = problem->nb_generated_col;
        CCutil_stop_timer(&(problem->tot_lb_lp_root), 0);

        switch (parms->bb_branch_strategy) {
            case conflict_strategy:
            case ahv_strategy:
                val = insert_into_branching_heap(root_pd, problem);
                CCcheck_val_2(val, "insert_into_branching_heap failed");
                break;

            case cbfs_conflict_strategy:
            case cbfs_ahv_strategy:
                insert_node_for_exploration(root_pd, problem);
                break;
        }
    }

    printf("GUB = %d, GLB = %d\n", problem->global_upper_bound,
           problem->global_lower_bound);

    if (problem->global_lower_bound != problem->global_upper_bound) {
        CCutil_start_resume_time(&(problem->tot_branch_and_bound));

        switch (parms->bb_branch_strategy) {
            case conflict_strategy:
                val = sequential_branching_conflict(problem);
                CCcheck_val(val, "Failed in sequential_branching_conflict");
                break;

            case cbfs_conflict_strategy:
                val = sequential_cbfs_branch_and_bound_conflict(problem);
                CCcheck_val_2(val, "Failed in CBFS conflict branching");
                break;
        }

        CCutil_stop_timer(&(problem->tot_branch_and_bound), 0);
        printf("Compute schedule finished with LB %d and UB %d\n",
               root_pd->lower_bound, problem->global_upper_bound);
    } else {
        problem->found = 1;
    }

    if (root_pd->lower_bound == problem->global_upper_bound) {
        problem->global_lower_bound = root_pd->lower_bound;
        problem->rel_error = (double)(problem->global_upper_bound -
                                      problem->global_lower_bound) /
                             ((double)problem->global_lower_bound);
        problem->status = optimal;
        printf("The optimal schedule is given by:\n");
        print_schedule(root_pd->bestcolors, root_pd->nbbest);
        printf("with total weighted completion time %d\n",
               root_pd->upper_bound);
    } else {
        problem->global_lower_bound = root_pd->lower_bound;
        problem->rel_error = (double)(problem->global_upper_bound -
                                      problem->global_lower_bound) /
                             ((double)problem->global_lower_bound);
        problem->status = meta_heur;
        problem->global_lower_bound = root_pd->lower_bound;
        printf("The suboptimal schedule is given by:\n");
        print_schedule(root_pd->bestcolors, root_pd->nbbest);
        printf("with total weighted completion time\n");
    }

    children_data_free(&problem->root_pd);
CLEAN:
    return val;
}
