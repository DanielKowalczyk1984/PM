#include <localsearch.h>
#include <wct.h>

static int add_feasible_solution(wctproblem *problem, solution *new_sol);
static int solution_set_c(solution *sol);
static void perturb_swap(solution *         sol,
                         local_search_data *data,
                         int                l1,
                         int                l2,
                         GRand *            rand_uniform);
void Perturb(solution *sol, local_search_data *data, GRand *rand_uniform);
void permutation_solution(GRand *rand_uniform, solution *sol);

int _job_compare_spt(const void *a, const void *b);
int compare_completion_time(BinomialHeapValue a, BinomialHeapValue b);
int compare_nb_job(gconstpointer a, gconstpointer b);

/**
 * comparefunctions
 */

int compare_completion_time(BinomialHeapValue a, BinomialHeapValue b) {
    partlist *x = (partlist *)a;
    partlist *y = (partlist *)b;
    int       C_a = x->c;
    int       C_b = y->c;
    int       key_a = x->key;
    int       key_b = y->key;

    if (C_a < C_b) {
        return -1;
    } else if (C_a > C_b) {
        return 1;
    } else if (key_a < key_b) {
        return -1;
    } else {
        return 1;
    }
}

int compare_nb_job(gconstpointer a, gconstpointer b){
    const Job *x = * (Job * const *) a;
    const Job *y = * (Job *const *) b;

    return (x->job - y ->job);
}

int _job_compare_spt(const void *a, const void *b) {
    const Job *x = ((const Job *)&a);
    const Job *y = ((const Job *)&b);

    if (x->processingime > y->processingime) {
        return 1;
    } else if (x->processingime < y->processingime) {
        return -1;
    } else if (x->duetime > y->duetime) {
        return 1;
    } else if (x->duetime < y->duetime) {
        return -1;
    } else if (x->weight > y->weight) {
        return 1;
    } else if (x->weight < y->weight) {
        return -1;
    } else if (x->job > y->job) {
        return 1;
    } else if (x->job < y->job) {
        return -1;
    }

    return (0);
}

/**
 * greedy constructions
 */

static int solution_set_c(solution *sol) {
    int           val = 0;
    partlist *    tmp = (partlist *)NULL;
    Job *         j = (Job *)NULL;
    BinomialHeap *heap =
        binomial_heap_new(BINOMIAL_HEAP_TYPE_MIN, compare_completion_time);
    CCcheck_NULL_2(heap, "Failed to allocate memory to heap");
    sol->tw = 0;
    sol->b = 0;

    for (int i = 0; i < sol->nmachines; ++i) {
        sol->part[i].c = 0;
        sol->part[i].tw = 0;
        sol->part[i].key = i;
        g_ptr_array_free(sol->part[i].machine, TRUE);
        sol->part[i].machine = g_ptr_array_new();
        binomial_heap_insert(heap, sol->part + i);
    }

    for (int i = 0; i < sol->njobs; ++i) {
        j = sol->perm[i];
        tmp = (partlist *)binomial_heap_pop(heap);
        j->index = tmp->machine->len;
        g_ptr_array_add(tmp->machine, j);
        tmp->c += j->processingime;
        sol->c[j->job] = tmp->c;
        tmp->tw += value_Fj(tmp->c, j);
        sol->tw += value_Fj(tmp->c, j);
        sol->b += j->duetime * (sol->njobs - i);
        binomial_heap_insert(heap, tmp);
    }

CLEAN:

    if (val) {
        solution_free(&sol);
    }

    binomial_heap_free(heap);
    return val;
}

int construct_spt(wctproblem *prob, solution *sol) {
    int val = 0;

    g_ptr_array_foreach(prob->g_job_array, g_set_sol_perm, sol);

    sol->njobs = prob->njobs;
    sol->nmachines = prob->nmachines;
    qsort(sol->perm, sol->njobs, sizeof(Job *), _job_compare_spt);
    val = solution_set_c(sol);
    CCcheck_val_2(val, "Failed in solution_set_c");
CLEAN:
    return val;
}

int construct_edd(wctproblem *prob, solution *sol) {
    int val = 0;

    g_ptr_array_foreach(prob->g_job_array, g_set_sol_perm, sol);

    sol->njobs = prob->njobs;
    sol->nmachines = prob->nmachines;
    val = solution_set_c(sol);
    CCcheck_val_2(val, "failed in solution_set_c");
CLEAN:
    return val;
}

int construct_random(wctproblem *prob, solution *sol, GRand *rand_uniform) {
    int val = 0;

    g_ptr_array_foreach(prob->g_job_array, g_set_sol_perm, sol);
    sol->njobs = prob->njobs;
    sol->nmachines = prob->nmachines;
    permutation_solution(rand_uniform, sol);
    val = solution_set_c(sol);
    CCcheck_val_2(val, "failed in solution_set_c");
CLEAN:
    return val;
}

void permutation_solution(GRand *rand_uniform, solution *sol) {
    int  i;
    Job *tmp = (Job *)NULL;

    for (i = 0; i <= sol->njobs - 2; i++) {
        int j = g_rand_int_range(rand_uniform, 0, sol->njobs - i);
        CC_SWAP(sol->perm[i], sol->perm[i + j], tmp);
    }
}

void RVND(solution *sol, local_search_data *data) {
    do {
        local_search_forward_insertion(sol, data, 1);

        if (data->updated) {
            continue;
        }

        local_search_swap_intra(sol, data, 0, 1);

        if (data->updated) {
            continue;
        }

        local_search_forward_insertion(sol, data, 2);

        if (data->updated) {
            continue;
        }

        local_search_swap_intra(sol, data, 0, 2);

        if (data->updated) {
            continue;
        }

        local_search_swap_intra(sol, data, 1, 1);

        if (data->updated) {
            continue;
        }

        local_search_insertion_inter(sol, data, 1);

        if (data->updated) {
            continue;
        }

        local_search_insertion_inter(sol, data, 2);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 1, 1);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 1, 2);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 1, 3);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 2, 2);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 2, 3);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 3, 3);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 3, 4);

        if (data->updated) {
            continue;
        }

        local_search_swap_inter(sol, data, 4, 4);

        if (data->updated) {
            continue;
        }
    } while (data->updated);
}

static void perturb_swap(solution *         sol,
                         local_search_data *data,
                         int                l1,
                         int                l2,
                         GRand *            rand_uniform) {
    int       m1, m2;
    unsigned  i1 = 0, i2 = 0;
    int       nmachines = sol->nmachines;
    Job **    tmp1 = (Job **)NULL;
    Job **    tmp2 = (Job **)NULL;
    Job *     tmp;
    partlist *part1 = (partlist *)NULL;
    partlist *part2 = (partlist *)NULL;
    m1 = g_rand_int_range(rand_uniform, 0, nmachines);
    m2 = g_rand_int_range(rand_uniform, 0, nmachines);

    while (m1 == m2) {
        m2 = g_rand_int_range(rand_uniform, 0, nmachines);
    }

    part1 = sol->part + m1;
    part2 = sol->part + m2;

    if (part1->machine->len <= (guint)l1 || part2->machine->len <= (guint)l2) {
        return;
    }

    tmp1 = CC_SAFE_MALLOC(l1, Job *);
    tmp2 = CC_SAFE_MALLOC(l2, Job *);
    sol->tw -= part1->tw + part2->tw;
    part1->c = 0;
    part1->tw = 0;
    part2->c = 0;
    part2->tw = 0;

    if (part1->machine->len - l1 != 0) {
        i1 = g_random_int_range(0, part1->machine->len - l1);
        assert(i1 < part1->machine->len - l1);
    }

    for (int i = 0; i < l1; ++i) {
        tmp1[i] = (Job *)g_ptr_array_index(part1->machine, i1);
        g_ptr_array_remove_index(part1->machine, i1);
    }

    if (part2->machine->len - l2 != 0) {
        i2 = (int)g_random_int_range(0, part2->machine->len - l2);
        assert(i2 < part2->machine->len - l2);
    }

    for (int i = 0; i < l2; ++i) {
        tmp2[i] = (Job *)g_ptr_array_index(part2->machine, i2);
        g_ptr_array_remove_index(part2->machine, i2);
    }

    if (part2->machine->len != 0) {
        i2 = g_random_int_range(0, part2->machine->len);
        assert(i2 < part2->machine->len);

        for (int i = 0; i < l1; ++i) {
            g_ptr_array_insert(part2->machine, i2 + i, tmp1[i]);
        }
    } else {
        for (int i = 0; i < l1; ++i) {
            g_ptr_array_add(part2->machine, tmp1[i]);
        }
    }

    if (part1->machine->len != 0) {
        i1 = g_random_int_range(0, part1->machine->len);
        assert(i1 < part1->machine->len);

        for (int i = 0; i < l2; ++i) {
            g_ptr_array_insert(part1->machine, i1 + i, tmp2[i]);
        }
    } else {
        for (int i = 0; i < l2; ++i) {
            g_ptr_array_add(part1->machine, tmp2[i]);
        }
    }

    for (unsigned i = 0; i < part1->machine->len; ++i) {
        tmp = (Job *)g_ptr_array_index(part1->machine, i);
        tmp->index = i;
        part1->c += tmp->processingime;
        sol->c[tmp->job] = part1->c;
        part1->tw += tmp->weight * CC_MAX(0, sol->c[tmp->job] - tmp->duetime);
    }

    for (unsigned i = 0; i < part2->machine->len; ++i) {
        tmp = (Job *)g_ptr_array_index(part2->machine, i);
        tmp->index = i;
        part2->c += tmp->processingime;
        sol->c[tmp->job] = part2->c;
        part2->tw += tmp->weight * CC_MAX(0, sol->c[tmp->job] - tmp->duetime);
    }

    sol->tw += part1->tw + part2->tw;
    CC_IFFREE(tmp1, Job *);
    CC_IFFREE(tmp2, Job *);
}

void Perturb(solution *sol, local_search_data *data, GRand *rand_uniform) {
    int L;
    L = g_rand_int_range(rand_uniform, 0, 3);

    for (int i = 0; i < L; ++i) {
        perturb_swap(sol, data, 1, 2, rand_uniform);
    }

    L = g_rand_int_range(rand_uniform, 0, 3);

    for (int i = 0; i < L; ++i) {
        perturb_swap(sol, data, 1, 3, rand_uniform);
    }

    L = g_rand_int_range(rand_uniform, 0, 3);

    for (int i = 0; i < L; ++i) {
        perturb_swap(sol, data, 2, 2, rand_uniform);
    }

    L = g_rand_int_range(rand_uniform, 0, 3);

    for (int i = 0; i < L; ++i) {
        perturb_swap(sol, data, 2, 3, rand_uniform);
    }

    for (int i = 0; i < sol->nmachines; ++i) {
        sol->part[i].used = 1;
    }

    local_search_create_W(sol, data);
    local_search_create_g(sol, data);
}

int heuristic_rpup(wctproblem *prob) {
    int    val = 0;
    GRand *rand_uniform = g_rand_new_with_seed(2011);
    g_random_set_seed(1984);
    int          ILS = prob->njobs;
    int          IR = 1;
    solution *   sol;
    solution *   sol1 = (solution *)NULL;
    CCutil_timer test;
    CCutil_init_timer(&test, (char *)NULL);
    local_search_data *data = (local_search_data *)NULL;
    local_search_data *data_RS = (local_search_data *)NULL;

    sol = solution_alloc(prob->nmachines, prob->njobs, prob->off);
    CCcheck_NULL_2(sol, "Failed to allocate memory");
    val = construct_edd(prob, sol);
    solution_print(sol);
    solution_canonical_order(sol, prob->e);
    CCcheck_val_2(val, "Failed construct edd");
    data = local_search_data_init(sol);
    CCcheck_NULL_2(data, "Failed to allocate memory to data");
    local_search_create_W(sol, data);
    local_search_create_g(sol, data);
    RVND(sol, data);
    prob->opt_sol = solution_alloc(prob->nmachines, prob->njobs, prob->off);
    CCcheck_NULL_2(prob->opt_sol, "Failed to allocate memory");
    solution_update(prob->opt_sol, sol);

    for (int i = 0; i < IR; ++i) {
        sol1 = solution_alloc(prob->nmachines, prob->njobs, prob->off);
        CCcheck_NULL_2(sol1, "Failed to allocate memory");
        val = construct_random(prob, sol1, rand_uniform);
        CCcheck_val_2(val, "Failed in construct random solution");
        data_RS = local_search_data_init(sol1);
        local_search_create_W(sol1, data_RS);
        local_search_create_g(sol1, data_RS);
        solution_update(sol, sol1);

        for (int j = 0; j < ILS; ++j) {
            RVND(sol1, data_RS);

            if (sol1->tw < sol->tw) {
                solution_update(sol, sol1);
                j = 0;
            }

            solution_update(sol1, sol);
            Perturb(sol1, data_RS, rand_uniform);
        }

        if (sol->tw < prob->opt_sol->tw) {
            solution_update(prob->opt_sol, sol);
        }

        local_search_data_free(&data_RS);
        solution_free(&sol1);
    }

    solution_print(prob->opt_sol);
    solution_canonical_order(prob->opt_sol, prob->e);
    for(int i = 0; i < prob->opt_sol->nmachines; ++i) {
        GPtrArray *machine = prob->opt_sol->part[i].machine;
        for(unsigned j = 0; j < machine->len; ++j) {
            Job *tmp = (Job *) g_ptr_array_index(machine, j);
            printf("%d ", prob->opt_sol->u[tmp->job]);
        }
        printf("\n");
    }
    solution_print(prob->opt_sol);
CLEAN:
    solution_free(&sol);
    local_search_data_free(&data);
    g_rand_free(rand_uniform);
    return val;
}

/** Construct feasible solutions */

void update_bestschedule(wctproblem *problem, solution *new_sol) {
    if (new_sol == NULL) {
        return;
    }

    if (new_sol->tw < problem->global_upper_bound) {
        problem->global_upper_bound = new_sol->tw;
        problem->rel_error = (double)(problem->global_upper_bound -
                                      problem->global_lower_bound) /
                             (problem->global_lower_bound);
        partlist_to_scheduleset(new_sol->part, new_sol->nmachines,
                                new_sol->njobs, &(problem->bestschedule),
                                &(problem->nbestschedule));
    }

    if (problem->global_upper_bound == problem->global_lower_bound) {
        problem->status = optimal;
    }
}

static int add_feasible_solution(wctproblem *problem, solution *new_sol) {
    int      val = 0;
    wctdata *root_pd = &(problem->root_pd);
    update_bestschedule(problem, new_sol);

    if (root_pd->ccount == 0 && problem->parms.construct != 0) {
        update_schedulesets(&root_pd->cclasses, &root_pd->ccount,
                            problem->bestschedule, problem->nbestschedule);
        root_pd->gallocated = root_pd->ccount;
    } else if (problem->parms.construct != 0) {
        partlist_to_scheduleset(new_sol->part, new_sol->nmachines,
                                new_sol->njobs, &(root_pd->newsets),
                                &(root_pd->nnewsets));
        add_newsets(root_pd);
    }

    return val;
}

int construct_feasible_solutions(wctproblem *problem) {
    int           val = 0;
    int           iterations = 0;
    wctdata *     pd = &(problem->root_pd);
    wctparms *    parms = &(problem->parms);
    CCutil_timer *timer = &(problem->tot_scatter_search);
    GRand *       rand1 = g_rand_new_with_seed(1984);
    GRand *       rand2 = g_rand_new_with_seed(1654651);
    CCutil_start_timer(timer);
    // while (1) {
    iterations++;
    solution *new_sol = solution_alloc(pd->nmachines, pd->njobs, problem->off);
    CCcheck_NULL(new_sol, "Failed to allocate")

        if (problem->status == no_sol) {}
    else {
        if (g_rand_boolean(rand1)) {
        } else {
        }
    }

    val = add_feasible_solution(problem, new_sol);
    CCcheck_val(val, "Failed in add_feasible_solution");
    // break;
    //}
    CCutil_suspend_timer(timer);
    printf("We needed %f seconds to construct %d solutions in %d iterations\n",
           timer->cum_zeit, parms->nb_feas_sol, iterations);
    printf("upperbound = %d, lowerbound = %d\n", problem->global_upper_bound,
           problem->global_lower_bound);
    CCutil_resume_timer(timer);
    g_rand_free(rand1);
    g_rand_free(rand2);
    CCutil_stop_timer(&(problem->tot_scatter_search), 0);
    return val;
}
