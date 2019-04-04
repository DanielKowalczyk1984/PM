#include <wct.h>

int debug = 0;

/*Information about debug*/
int  dbg_lvl() { return debug; }
void set_dbg_lvl(int dbglvl) { debug = dbglvl; }

/*Functions for initialization of the problem and freeing the problem*/
void wctproblem_init(Problem *problem) {
    /** Job data */
    problem->g_job_array = g_ptr_array_new_with_free_func(g_job_free);
    problem->opt_sol = (Solution *)NULL;
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
    /*B&B info*/
    problem->nwctdata = 0;
    problem->global_upper_bound = INT_MAX;
    problem->root_lower_bound = 0;
    problem->global_lower_bound = 0;
    problem->root_upper_bound = INT_MAX;
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
    wctdata_init(&(problem->root_pd), problem);
    set_id_and_name(&(problem->root_pd), 0, "root_node");
    problem->nwctdata++;
    /*parms of the problem*/
    parms_init(&(problem->parms));
    /*heap initialization*/
    problem->unexplored_states = g_ptr_array_new();
    problem->non_empty_level_pqs = g_queue_new();
    problem->last_explored = -1;
    problem->found = 0;
    problem->nb_explored_nodes = 0;
    problem->nb_generated_col = 0;
    /*CPU timer initialisation*/
    CCutil_init_timer(&(problem->tot_build_dd), "tot_build_dd");
    CCutil_init_timer(&(problem->tot_cputime), "tot_cputime");
    CCutil_init_timer(&(problem->tot_branch_and_bound), "tot_branch_and_bound");
    CCutil_init_timer(&(problem->tot_strong_branching), "tot_strong_branching");
    CCutil_init_timer(&(problem->tot_lb_root), "tot_lb_root");
    CCutil_init_timer(&(problem->tot_lb), "tot_lb");
    CCutil_init_timer(&(problem->tot_solve_lp), "tot_solve_lp");
    CCutil_init_timer(&(problem->tot_pricing), "tot_pricing");
    CCutil_init_timer(&(problem->tot_heuristic), "tot_heuristic");
    /** initialize colPool */
    problem->ColPool = g_ptr_array_new_with_free_func(g_scheduleset_free);
    problem->nArtificials = 0;
    /** initialize the time */
    problem->real_time_build_dd = 0.0;
    problem->real_time_total = getRealTime();
    problem->real_time_branch_and_bound = 0.0;
    problem->real_time_strong_branching = 0.0;
    problem->real_time_lb_root = 0.0;
    problem->real_time_lb = 0.0;
    problem->real_time_pricing = 0.0;
    problem->real_time_heuristic = 0.0;
    CCutil_start_timer(&(problem->tot_cputime));
}

void wctproblem_free(Problem *problem) {
    /*free the parameters*/
    parms_free(&(problem->parms));
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
    g_ptr_array_free(problem->unexplored_states, TRUE);
    g_ptr_array_free(problem->ColPool, TRUE);
    g_queue_free(problem->non_empty_level_pqs);
    schedulesets_free(&(problem->initsets), &(problem->gallocated));
    schedulesets_free(&(problem->bestschedule), &(problem->nbestschedule));
    solution_free(&(problem->opt_sol));
}

/*Functions for initialization and free the data*/
void wctdata_init(wctdata *pd, Problem *prob) {
    /*Initialization B&B data*/
    pd->id = -1;
    pd->depth = 0;
    pd->status = initialized;
    sprintf(pd->pname, "temporary");
    /*Initialization node instance data*/
    pd->njobs = 0;
    pd->orig_node_ids = (int *)NULL;
    pd->H_max = 0;
    pd->H_min = 0;
    pd->local_intervals = g_ptr_array_new_with_free_func(g_interval_free);
    pd->ordered_jobs = g_ptr_array_new_with_free_func(g_free);
    pd->jobarray = (GPtrArray *) NULL;
    pd->sump = (int **)NULL;
    /** Initialization data */
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
    pd->MIP = (wctlp *)NULL;
    pd->x = (double *)NULL;
    pd->x_e = (double *) NULL;
    pd->coef = (double *)NULL;
    pd->pi = (double *)NULL;
    pd->kpc_pi = (int *)NULL;
    /**init stab data */
    pd->pi_in = (double *)NULL;
    pd->pi_out = (double *)NULL;
    pd->pi_sep = (double *)NULL;
    pd->subgradient = (double *)NULL;
    pd->subgradient_in = (double *)NULL;
    pd->reduced_cost = 0.0;
    pd->alpha = 0.8;
    pd->update = 1;
    pd->iterations = 0;
    pd->hasstabcenter = 0;
    pd->beta = 0.0;
    pd->node_stab = -1;
    pd->subgradientnorm = 0.0;
    pd->dualdiffnorm = 0.0;
    pd->hybridfactor = 0.0;
    pd->subgradientproduct = 0.0;
    /*Initialization pricing_problem*/
    pd->solver = (PricerSolver *)NULL;
    pd->nnonimprovements = 0;
    pd->dzcount = 0;
    pd->bestcolors = (scheduleset *)NULL;
    pd->nbbest = 0;
    pd->debugcolors = (scheduleset *)NULL;
    pd->ndebugcolors = 0;
    pd->opt_track = 0;
    pd->localColPool = g_ptr_array_new_with_free_func(g_scheduleset_free);
    pd->cstat = (int *)NULL;
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
    pd->v1 = (Job *)NULL;
    pd->v2 = (Job *)NULL;
    /** Wide branching */
    pd->v1_wide = (int *)NULL;
    pd->v2_wide = (int *)NULL;
    pd->nb_wide = 0;
    pd->same_children_wide = (wctdata **)NULL;
    pd->diff_children_wide = (wctdata **)NULL;

    pd->problem = prob;
}

void lpwctdata_free(wctdata *pd) {
    /**
     * free all the gurobi data associated with the LP relaxation
     */
    if (pd->LP) {
        wctlp_free(&(pd->LP));
    }

    /**
     * free all the data associated with the LP
     */
    CC_IFFREE(pd->coef, double);
    CC_IFFREE(pd->pi, double);
    CC_IFFREE(pd->x, double);
    CC_IFFREE(pd->x_e, double);
    CC_IFFREE(pd->kpc_pi, int);
    CC_IFFREE(pd->pi_out, double);
    CC_IFFREE(pd->pi_in, double);
    CC_IFFREE(pd->pi_sep, double);
    CC_IFFREE(pd->subgradient, double);
    CC_IFFREE(pd->subgradient_in, double);
    CC_IFFREE(pd->rhs, double);
    CC_IFFREE(pd->cstat, int);

    /**
     * free all the schedules from the localColPool
     */
    g_ptr_array_free(pd->localColPool,TRUE);
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
    // g_ptr_array_free(pd->localColPool, TRUE);
    if (pd->solver) {
        freeSolver(pd->solver);
        pd->solver = (PricerSolver *)NULL;
    }
}

void wctdata_free(wctdata *pd) {
    schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    temporary_data_free(pd);
    if (pd->sump) {
        for (unsigned i = 0; i < pd->local_intervals->len; ++i) {
            CC_IFFREE(pd->sump[i], int);
        }
        CC_IFFREE(pd->sump, int *)
    }

    g_ptr_array_free(pd->local_intervals, TRUE);
    g_ptr_array_free(pd->ordered_jobs, TRUE);
    CC_IFFREE(pd->elist_same, int);
    CC_IFFREE(pd->elist_differ, int);
    CC_IFFREE(pd->v1_wide, int);
    CC_IFFREE(pd->v2_wide, int);
    CC_IFFREE(pd->orig_node_ids, int);
    CC_IFFREE(pd->v1_wide, int);
    CC_IFFREE(pd->v2_wide, int);
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

static int prefill_heap(wctdata *pd, Problem *problem) {
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

int compute_schedule(Problem *problem) {
    int       val = 0;
    wctdata * root_pd = &(problem->root_pd);
    Parms *parms = &(problem->parms);
    problem->mult_key = 1.0;
    problem->root_upper_bound = problem->global_upper_bound;
    problem->root_lower_bound = problem->global_lower_bound;
    problem->root_rel_error =
        (double)(problem->global_upper_bound - problem->global_lower_bound) /
        ((double)problem->global_lower_bound);
    prune_duplicated_sets(root_pd);
    init_BB_tree(problem);
    print_size_to_csv(problem, root_pd);

    if (root_pd->status >= LP_bound_computed) {
        val = prefill_heap(root_pd, problem);
        CCcheck_val(val, "Failed in prefill_heap");
    } else {
        CCutil_start_timer(&(problem->tot_lb_root));
        val = compute_lower_bound(problem, root_pd);
        CCcheck_val_2(val, "Failed in compute_lower_bound");

        if (root_pd->lower_bound > problem->global_lower_bound) {
            problem->global_lower_bound = root_pd->lower_bound;
            problem->root_lower_bound = root_pd->lower_bound;
            problem->root_rel_error = (double)(problem->root_upper_bound -
                                               problem->root_lower_bound) /
                                      ((double)problem->root_lower_bound);
        }

        CCcheck_val_2(val, "Failed in compute_lower_bound");
        problem->nb_generated_col_root = problem->nb_generated_col;
        CCutil_stop_timer(&(problem->tot_lb_root), 0);

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

int add_solution_to_colpool(Solution *sol, wctdata *pd) {
    int          val = 0;

    for (int i = 0; i < sol->nmachines; ++i) {
        GPtrArray *machine = sol->part[i].machine;
        scheduleset *tmp = scheduleset_from_solution(machine, pd->njobs);
        CCcheck_NULL_2(tmp, "Failed to allocate memory");
        g_ptr_array_add(pd->localColPool, tmp);
    }

CLEAN:
    return val;
}

int add_solution_to_colpool_and_lp(Solution *sol, wctdata *pd) {
    int          val = 0;
    scheduleset *tmp;

    for (int i = 0; i < sol->nmachines; ++i) {
        GPtrArray *machine = sol->part[i].machine;
        tmp = scheduleset_from_solution(machine, pd->njobs);
        CCcheck_NULL_2(tmp, "Failed to allocate memory");
        g_ptr_array_add(pd->localColPool, tmp);
    }

    for (unsigned i = 0; i < pd->localColPool->len; ++i) {
        tmp = (scheduleset *)g_ptr_array_index(pd->localColPool, i);
        addColToLP(tmp, pd);
    }

CLEAN:
    return val;
}
