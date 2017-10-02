#include <wct.h>

static int collect_same_child_conflict(wctdata *cd);
static int collect_diff_child_conflict(wctdata *cd);
static int remove_finished_subtree_conflict(wctdata *child);
static int compare_nodes_bfs(BinomialHeapValue a, BinomialHeapValue b);
static int compare_nodes_dfs(BinomialHeapValue a, BinomialHeapValue b);
static int get_int_heap_key(double dbl_heap_key,
                            int    v1,
                            int    v2,
                            int    nmachines,
                            int    njobs,
                            double error);
static int get_int_heap_key_0(double dbl_heap_key, int v1, int v2);

static int compare_nodes_dfs(BinomialHeapValue a, BinomialHeapValue b) {
    wctdata *x = (wctdata *)a;
    wctdata *y = (wctdata *)b;
    double lp_a = (x->LP_lower_bound_BB - x->depth * 10000 - 0.5 * (x->id % 2));
    double lp_b = (y->LP_lower_bound_BB - y->depth * 10000 - 0.5 * (y->id % 2));

    if (lp_a < lp_b) {
        return -1;
    } else {
        return 1;
    }
}

static int compare_nodes_bfs(BinomialHeapValue a, BinomialHeapValue b) {
    double *lp_a = &(((wctdata *)a)->eta_in);
    double *lp_b = &(((wctdata *)b)->eta_in);

    if (*lp_a < *lp_b) {
        return -1;
    } else {
        return 1;
    }
}

MAYBE_UNUSED
static int x_frac(const double x, double error) {
    double mean = error;
    double frac = fabs(x - mean);
    assert(frac <= 1.0);
    return (int)(frac * (double)INT_MAX);
}

MAYBE_UNUSED
static int get_int_heap_key(double dbl_heap_key,
                            int    v1,
                            int    v2,
                            int    nmachines,
                            int    njobs,
                            double error) {
    int    val = INT_MAX - 1;
    double temp;

    if (dbl_heap_key > 0.5 && nmachines == 1) {
        temp = 1 - dbl_heap_key;
    } else {
        temp = dbl_heap_key;
    }

    double error2 = dbl_heap_key > error
                        ? (ABS(error - temp) + 0.0001) / (error)
                        : (ABS(error - temp) + 0.0001) / (error);
    error2 = error2 / njobs;

    if (dbl_heap_key >= error) {
        if (dbl_heap_key >= 1.0 - lp_int_tolerance()) {
            return val;
        }

        val = x_frac(MIN(1.0, temp + ABS((v2 - v1)) * error2), error);
    } else {
        val = x_frac(MAX(0.0, temp - ABS((v2 - v1)) * error2), error);
    }

    return val;
}

MAYBE_UNUSED
static int get_int_heap_key_0(double dbl_heap_key, int v1, int v2) {
    int val = INT_MAX - 1;

    if (dbl_heap_key >= 0.5) {
        if (dbl_heap_key >= 1.0) {
            return val;
        }

        return x_frac(dbl_heap_key / ABS((v2 - v1)), 0.5);
    }

    return x_frac(dbl_heap_key / ABS((v2 - v1)), 0.5);
}

void init_BB_tree(wctproblem *problem) {
    switch (problem->parms.bb_branch_strategy) {
        case conflict_strategy:
            problem->br_heap_a =
                binomial_heap_new(BINOMIAL_HEAP_TYPE_MIN, compare_nodes_dfs);
            break;

        case cbfs_conflict_strategy:
        case cbfs_ahv_strategy:
            problem->br_heap_a = (BinomialHeap *)NULL;
            break;
    }
}

int insert_frac_pairs_into_heap(wctdata *    pd,
                                int *        nodepair_refs,
                                double *     nodepair_weights,
                                int          npairs,
                                pmcheap *    heap) {
    int     val = 0;
    int     i;
    int     ref_key;
    int nb_cols;
    double *mean_error = CC_SAFE_MALLOC(npairs, double);
    int *   mean_counter = CC_SAFE_MALLOC(npairs, int);
    scheduleset *tmp_schedule;
    Job *tmp_j1;
    Job *tmp_j2;

    CCcheck_NULL_2(mean_error, "Failed to allocate memory");
    CCcheck_NULL_2(mean_counter, "Failed to allocate memory");
    fill_dbl(mean_error, npairs, 0.0);
    fill_int(mean_counter, npairs, 0);
    wctlp_get_nb_cols(pd->LP, &nb_cols);
    assert(nb_cols == pd->localColPool->len);


    for (i = 0; i < nb_cols; ++i) {
        int j;

        if (pd->x[i] <= 0.0 + lp_int_tolerance() ||
            pd->x[i] >= 1.0 - lp_int_tolerance()) {
            continue;
        }
        tmp_schedule = (scheduleset *) g_ptr_array_index(pd->localColPool, i);

        for(j = 0; j < tmp_schedule->partial_machine->len; ++j) {
            tmp_j1 = (Job *) g_ptr_array_index(tmp_schedule->partial_machine,j);
            int v1 = tmp_j1->job;
            int k;
            ref_key = nodepair_ref_key(v1, v1);
            nodepair_weights[ref_key] += pd->x[i];

            for (k = j + 1; k < tmp_schedule->partial_machine->len; ++k) {
                assert(k != j);
                tmp_j2 = (Job *) g_ptr_array_index(tmp_schedule->partial_machine,k);
                int v2 = tmp_j2->job;
                assert(v1 < v2);
                ref_key = nodepair_ref_key(v1, v2);
                mean_error[ref_key] += 0.5 - ABS(pd->x[i] - floor(pd->x[i]) - 0.5);
                mean_counter[ref_key]++;
                nodepair_weights[ref_key] += pd->x[i];
            }
        }
    }

    for (ref_key = 0; ref_key < npairs; ++ref_key) {
        int v1, v2;
        inodepair_ref_key(&v1, &v2, ref_key);

        if (v1 != v2 && nodepair_weights[ref_key] > 0.0) {
            int v1_key = nodepair_ref_key(v1, v1);
            int v2_key = nodepair_ref_key(v2, v2);
            mean_error[ref_key] = mean_error[ref_key] / mean_counter[ref_key];
            double denom =
                (nodepair_weights[v1_key] + nodepair_weights[v2_key]) / 2;
            double dbl_heap_key = nodepair_weights[ref_key] / denom;
            // printf("error %f %d %d %f\n",dbl_heap_key, v1, v2,
            // mean_error[ref_key]);
            int int_heap_key =
                get_int_heap_key(dbl_heap_key, v1, v2, mean_counter[ref_key],
                                 pd->njobs, mean_error[ref_key]);
            val = pmcheap_insert(heap, int_heap_key + 1,
                                 (void *)&(nodepair_refs[ref_key]));
            CCcheck_val_2(val, "Failed in pmcheap_insert");
        }
    }

    if (dbg_lvl()) {
        printf("Size of frac heap is %d\n", pmcheap_size(heap));
    }

CLEAN:
    CC_IFFREE(mean_counter, int)
    CC_IFFREE(mean_error, double)
    return val;
}

static int trigger_lb_changes_conflict(wctdata *child) {
    int      val = 0;
    int      i;
    int      new_lower_bound = child->lower_bound;
    wctdata *pd = (wctdata *)child->parent;

    while (pd) {
        for (i = 0; i < pd->nsame; ++i) {
            if (pd->same_children[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->same_children[i].lower_bound;
            }
        }

        for (i = 0; i < pd->ndiff; ++i) {
            if (pd->diff_children[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->diff_children[i].lower_bound;
            }
        }

        if (new_lower_bound > pd->lower_bound) {
            if (!pd->parent) { /* i.e. pd == root_cd */
                time_t current_time;
                char   current_timestr[40] = "";
                (void)time(&current_time);
                strftime(current_timestr, 39, "%c", localtime(&current_time));
                printf("Lower bound increased from %d to %d (%s). \n",
                       pd->lower_bound, new_lower_bound, current_timestr);
            }

            pd->lower_bound = new_lower_bound;
            pd = pd->parent;
        } else {
            pd = (wctdata *)NULL;
        }
    }

    return val;
}

void adapt_global_upper_bound(wctproblem *problem, int new_upper_bound) {
    if (problem->global_upper_bound > new_upper_bound) {
        problem->global_upper_bound = new_upper_bound;
    }
}

static int collect_same_child_conflict(wctdata *cd) {
    int rval = 0;
    int c;

    for (c = 0; c < cd->nsame; ++c) {
        if (cd->same_children[c].nbbest &&
            (!cd->nbbest ||
             cd->same_children[c].besttotwct < cd->upper_bound)) {
            if (cd->nbbest) {
                schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct = cd->same_children[c].besttotwct;
            cd->nbbest = cd->same_children[c].nbbest;
            cd->same_children[c].nbbest = 0;
            cd->bestcolors = cd->same_children[c].bestcolors;
            cd->same_children[c].bestcolors = (scheduleset *)NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_diff_child_conflict(wctdata *cd) {
    int rval = 0;
    int c;

    for (c = 0; c < cd->ndiff; ++c) {
        if (cd->diff_children[c].nbbest &&
            (!cd->nbbest ||
             cd->diff_children[c].besttotwct < cd->upper_bound)) {
            if (cd->nbbest) {
                schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct = cd->diff_children[c].besttotwct;
            cd->nbbest = cd->diff_children[c].nbbest;
            cd->diff_children[c].nbbest = 0;
            cd->bestcolors = cd->diff_children[c].bestcolors;
            cd->diff_children[c].bestcolors = (scheduleset *)NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int remove_finished_subtree_conflict(wctdata *child) {
    int      val = 0;
    int      i;
    wctdata *cd = (wctdata *)child;
    int      all_same_finished = 1;
    int      all_diff_finished = 1;

    while (cd) {
        for (i = 0; i < cd->nsame; ++i) {
            if (cd->same_children[i].status < infeasible) {
                all_same_finished = 0;
                break;
            }
        }

        if (cd->nsame && all_same_finished) {
            val = collect_same_child_conflict(cd);
            CCcheck_val_2(val, "Failed in collect_same_children");

            for (i = 0; i < cd->nsame; ++i) {
                wctdata_free(cd->same_children + i);
            }

            free(cd->same_children);
            cd->same_children = (wctdata *)NULL;
            cd->nsame = 0;
        }

        for (i = 0; i < cd->ndiff; ++i) {
            if (cd->diff_children[i].status < infeasible) {
                all_diff_finished = 0;
                break;
            }
        }

        if (cd->ndiff && all_diff_finished) {
            val = collect_diff_child_conflict(cd);
            CCcheck_val_2(val, "Failed in collect_diff_children");

            for (i = 0; i < cd->ndiff; ++i) {
                wctdata_free(cd->diff_children + i);
            }

            free(cd->diff_children);
            cd->diff_children = (wctdata *)NULL;
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

static void scheduleset_unify(GPtrArray *array) {
    int         i;
    int it = 1;
    int first_del = -1;
    int last_del = -1;
    int nb_col = array->len;
    scheduleset *temp, *prev;
    g_ptr_array_sort(array, g_scheduleset_less);

    if (!(array->len)) {
        return;
    }

    prev = (scheduleset *) g_ptr_array_index(array, 0);
    /* Find first non-empty set */
    for (i = 1; i < nb_col; ++i) {
        temp = (scheduleset *) g_ptr_array_index(array, it);
        if (scheduleset_less(prev,temp)) {
            if (first_del != -1) {
                /** Delete recently found deletion range.*/
                g_ptr_array_remove_range(array, first_del, last_del - first_del + 1);
                it = it - (last_del - first_del);
                first_del = last_del = -1;
            } else{
                it++;
            }
            prev = temp;
        } else {
            if (first_del == -1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            prev = temp;
            it++;
        }
    }

    if (first_del != -1) {
        g_ptr_array_remove_range(array, first_del, last_del - first_del + 1);
    }
}

int prune_duplicated_sets(wctdata *pd) {
    int val = 0;
    scheduleset *tmp;
    GPtrArray *tmp_a;
    Job *tmp_j;
    scheduleset_unify(pd->localColPool);

    if (dbg_lvl() > 1) {
        for (int i = 0; i < pd->localColPool->len; ++i) {
            tmp = (scheduleset *) g_ptr_array_index(pd->localColPool, i);
            tmp_a = tmp->partial_machine;

            printf("TRANSSORT SET ");

            for (int j = 0; j < tmp_a->len; ++j) {
                tmp_j = (Job *) g_ptr_array_index(tmp_a, j);
                printf(" %d", tmp_j->job);
            }

            printf("\n");
        }
    }

    return val;
}

int skip_wctdata(wctdata *pd, wctproblem *problem) {
    BinomialHeap *br_heap = problem->br_heap_a;

    if (dbg_lvl() > 0) {
        printf(
            "Skipping with lb %d and ub %d at depth %d (id = %d, "
            "opt_track = %d, unprocessed nodes = %u).\n",
            pd->lower_bound, pd->upper_bound, pd->depth, pd->id, pd->opt_track,
            binomial_heap_num_entries(br_heap));
    }

    pd->status = finished;
    return 0;
}

int insert_into_branching_heap(wctdata *pd, wctproblem *problem) {
    int       val = 0;
    int       heap_key = 0;
    wctparms *parms = &(problem->parms);
    problem->parms.bb_search_strategy = dfs_strategy;
    int lb = pd->lower_bound;

    if (dbg_lvl()) {
        printf(
            "Inserting into branching heap with lb %d and ub %d at depth %d "
            "(id "
            "= %d) heap_key = %d\n",
            pd->lower_bound, pd->upper_bound, pd->depth, pd->id, heap_key);
    }

    free_elist(pd, &(problem->parms));

    if (lb < pd->upper_bound && pd->status != infeasible) {
        binomial_heap_insert(problem->br_heap_a, pd);
        CCcheck_val(val, "Failed at pmcheap_insert");
    } else {
        skip_wctdata(pd, problem);
    }

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
            val = trigger_lb_changes_conflict(pd);
            CCcheck_val_2(val, "Failed in trigger_lb_changes_conflict");
            break;
    }

CLEAN:
    return val;
}

void insert_node_for_exploration(wctdata *pd, wctproblem *problem) {
    unsigned int level = pd->ecount_same;
    wctparms *   parms = &(problem->parms);

    if (pd->lower_bound < pd->upper_bound && pd->status != infeasible) {
        while (problem->unexplored_states->len <= level) {
            g_ptr_array_add(problem->unexplored_states,
                            (gpointer)binomial_heap_new(BINOMIAL_HEAP_TYPE_MIN,
                                                        compare_nodes_bfs));
        }

        unsigned int wasPrevEmpty =
            binomial_heap_num_entries((
                BinomialHeap *)(problem->unexplored_states->pdata[level])) == 0;
        binomial_heap_insert(
            (BinomialHeap *)(problem->unexplored_states->pdata[level]), pd);

        if (wasPrevEmpty) {
            if (level == problem->last_explored) {
                g_queue_push_tail(problem->non_empty_level_pqs,
                                  (problem->unexplored_states->pdata[level]));
                return;
            }

            unsigned int isPrevLevelEmpty =
                (level > 0)
                    ? binomial_heap_num_entries(
                          (BinomialHeap *)
                              problem->unexplored_states->pdata[level - 1]) == 0
                    : TRUE;

            if (isPrevLevelEmpty) {
                g_queue_push_head(problem->non_empty_level_pqs,
                                  (problem->unexplored_states->pdata[level]));
            } else {
                gpointer data = g_queue_pop_head(problem->non_empty_level_pqs);
                g_queue_push_head(problem->non_empty_level_pqs,
                                  (problem->unexplored_states->pdata[level]));
                g_queue_push_head(problem->non_empty_level_pqs, data);
            }
        }
    } else {
        skip_wctdata(pd, problem);
    }

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
        case cbfs_conflict_strategy:
            trigger_lb_changes_conflict(pd);
            break;
    }
}

wctdata *get_next_node(wctproblem *problem) {
    wctdata *     pd = (wctdata *)NULL;
    GQueue *      non_empty_level_pqs = problem->non_empty_level_pqs;
    BinomialHeap *next_level_pq =
        (BinomialHeap *)g_queue_pop_head(non_empty_level_pqs);

    if (next_level_pq == (BinomialHeap *)NULL) {
        return pd;
    }

    pd = (wctdata *)binomial_heap_pop(next_level_pq);

    while (pd->lower_bound >= problem->global_upper_bound) {
        if (binomial_heap_num_entries(next_level_pq) == 0) {
            if (non_empty_level_pqs->length == 0) {
                problem->last_explored = pd->ecount_same;
                return pd;
            }

            next_level_pq =
                (BinomialHeap *)g_queue_pop_head(non_empty_level_pqs);
        }

        pd = (wctdata *)binomial_heap_pop(next_level_pq);
    }

    if (binomial_heap_num_entries(next_level_pq) > 0) {
        g_queue_push_tail(non_empty_level_pqs, next_level_pq);
    }

    problem->last_explored = pd->ecount_same;
    return pd;
}

void free_elist(wctdata *cd, wctparms *parms) {
    if (cd->parent && parms->delete_elists) {
        CC_IFFREE(cd->elist_same, int);
        CC_IFFREE(cd->elist_differ, int);
        CC_IFFREE(cd->v1_wide, int);
        CC_IFFREE(cd->v2_wide, int);
        CC_IFFREE(cd->orig_node_ids, int);
    }
}

int branching_msg(wctdata *pd, wctproblem *problem) {
    BinomialHeap *heap = problem->br_heap_a;
    wctdata *     root = &(problem->root_pd);

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf(
            "Branching with lb %d (LP %f) at depth %d (id = %d, "
            "time = %f, unprocessed nodes = %u, nbjobs= %d, upper bound = %d, "
            "lower bound = %d, v1 = %d, v2 = %d, nbdiff = %d, nbsame = %d, ZDD "
            "size= %zu, nb_cols = %u ).\n",
            pd->lower_bound, pd->LP_lower_bound_BB, pd->depth, pd->id,
            problem->tot_cputime.cum_zeit, binomial_heap_num_entries(heap),
            pd->njobs, problem->global_upper_bound, root->lower_bound, pd->v1->job,
            pd->v2->job, pd->ecount_differ, pd->ecount_same,
            get_datasize(pd->solver), pd->localColPool->len);
        CCutil_resume_timer(&problem->tot_cputime);
        problem->nb_explored_nodes++;
    }

    return 0;
}

int sequential_branching_conflict(wctproblem *problem) {
    int           val = 0;
    wctdata *     pd;
    BinomialHeap *br_heap = problem->br_heap_a;
    wctparms *    parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL BRANCHING CONFLICT:\n");
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
            remove_finished_subtree_conflict(pd);
        } else {
            branching_msg(pd, problem);

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches_conflict(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nsame; i++) {
                val = insert_into_branching_heap(&(pd->same_children[i]),
                                                 problem);
                CCcheck_val_2(val, "Failed in insert_into_branching_heap");
            }

            for (i = 0; i < pd->ndiff; i++) {
                val =
                    insert_into_branching_heap(pd->diff_children + i, problem);
                CCcheck_val_2(val, "Faield at insert_into_branching_heap");
            }

            // assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                if (pd->depth == 0) {
                    problem->found = 1;
                }

                remove_finished_subtree_conflict(pd);
            }
        }

        CCutil_suspend_timer(&problem->tot_branch_and_bound);
    }

    CCutil_resume_timer(&problem->tot_branch_and_bound);

    if (pd) {
        printf("Branching timeout of %f second reached\n",
               parms->branching_cpu_limit);
    }

CLEAN:
    return val;
}

int branching_msg_cbfs(wctdata *pd, wctproblem *problem) {
    int      nb_nodes = 0;
    wctdata *root = &(problem->root_pd);

    for (unsigned int i = 0; i < problem->unexplored_states->len; ++i) {
        BinomialHeap *heap =
            (BinomialHeap *)(problem->unexplored_states->pdata[i]);
        nb_nodes += binomial_heap_num_entries(heap);
    }

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf(
            "Branching with lb %d (LP %f) at depth %d (id = %d, "
            "time = %f, unprocessed nodes = %d, nbjobs= %d, upper bound = %d, "
            "lower bound = %d, v1 = %d, v2 = %d, nbdiff = %d, nbsame = %d, ZDD "
            "size = %zu, nb_cols = %u ).\n",
            pd->lower_bound, pd->LP_lower_bound, pd->depth, pd->id,
            problem->tot_cputime.cum_zeit, nb_nodes, pd->njobs,
            problem->global_upper_bound, root->lower_bound, pd->v1->job, pd->v2->job,
            pd->ecount_differ, pd->ecount_same, get_datasize(pd->solver),
            pd->localColPool->len);
        CCutil_resume_timer(&problem->tot_cputime);
        problem->nb_explored_nodes++;
    }

    return 0;
}

int sequential_cbfs_branch_and_bound_conflict(wctproblem *problem) {
    int       val = 0;
    wctdata * pd;
    wctparms *parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL BRANCHING CONFLICT + CBFS SEARCHING:\n");
    CCutil_suspend_timer(&problem->tot_branch_and_bound);

    while ((pd = get_next_node(problem)) &&
           problem->tot_branch_and_bound.cum_zeit <
               parms->branching_cpu_limit) {
        CCutil_resume_timer(&problem->tot_branch_and_bound);
        int i;
        pd->upper_bound = problem->global_upper_bound;

        if (pd->lower_bound >= pd->upper_bound ||
            pd->eta_in > pd->upper_bound - 1) {
            skip_wctdata(pd, problem);
            remove_finished_subtree_conflict(pd);
        } else {
            branching_msg_cbfs(pd, problem);

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches_conflict(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nsame; i++) {
                insert_node_for_exploration(pd->same_children + i, problem);
            }

            for (i = 0; i < pd->ndiff; i++) {
                insert_node_for_exploration(pd->diff_children + i, problem);
            }

            // assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                if (pd->depth == 0) {
                    problem->found = 1;
                }

                remove_finished_subtree_conflict(pd);
            }
        }

        CCutil_suspend_timer(&problem->tot_branch_and_bound);
    }

    CCutil_resume_timer(&problem->tot_branch_and_bound);

    if (pd) {
        printf("Branching timeout of %f second reached\n",
               parms->branching_cpu_limit);
    }

CLEAN:
    return val;
}
