#include "wct.h"

static int add_feasible_solution(wctproblem *problem, solution *new_sol);
static int solution_set_c(solution *sol);

/**
 * comparefunctions
 */

int compare_func1(gconstpointer a, gconstpointer b, void *user_data) {
    const int *v = &(((const partlist *)a)->c);
    const int *w = &(((const partlist *)b)->c);

    if (*v != *w) {
        return *v - *w;
    } else {
        if (*v == 0 || *w == 0) {
            return *v - *w;
        }

        const int *vv = &(((Job *)((const partlist *)a)->list->head->data)->job);
        const int *ww = &(((Job *)((const partlist *)b)->list->head->data)->job);
        return *vv - *ww;
    }
}

int compare_completion_time(BinomialHeapValue a, BinomialHeapValue b) {
    partlist *x = (partlist *)a;
    partlist *y = (partlist *)b;
    int C_a = x->c;
    int C_b = y->c;
    int key_a = x->key;
    int key_b = y->key;

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

int _job_compare_spt(const void *a, const void *b) {
    Job *x = *((Job **) a);
    Job *y = *((Job **) b);

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

int random_rcl_assignment(Job *jobarray, int njobs, int nmachines,
                          solution *new_sol, GRand *rand_) {
    int i;
    double max;
    double min;
    double lb, ub;
    double temp_dbl;
    partlist *temp = (partlist *) NULL;
    Job *temp_job = (Job *) NULL;
    //GList *it = (GList *) NULL;
    GQueue *to_do_list = (GQueue *) NULL;
    temp = new_sol->part;
    to_do_list = g_queue_new();

    for (i = 0; i < njobs; ++i) {
        g_queue_push_tail(to_do_list, jobarray + i);
    }

    while (!g_queue_is_empty(to_do_list)) {
        temp_job = (Job *)to_do_list->head->data;
        max = ((double)temp[0].c + (double)temp_job->processingime);
        min = max;
        GArray *rcl = g_array_new(FALSE, FALSE, sizeof(pair_job_machine));

        /** Compute min and max */
        for (i = 1; i < nmachines; ++i) {
            //for (it = to_do_list->head; it; it = it->next)
            //{
            //temp_job = (Job*)it->data;
            temp_dbl = (temp[i].c + temp_job->processingime);

            if (max < temp_dbl) {
                max = temp_dbl;
            }

            if (min > temp_dbl) {
                min = temp_dbl;
            }

            //}
        }

        /** Compute RCL */
        pair_job_machine temp_job_machine;
        lb = min;
        ub = min + 0.25 * (max - lb);

        for (i = 0; i < nmachines; ++i) {
            //for (it = to_do_list->head; it; it = g_list_next(it))
            //{
            //temp_job = ((Job*)it->data);
            double g = ((double)temp[i].c + (double)temp_job->processingime);

            if (lb <= g && g <= ub) {
                temp_job_machine.job = temp_job->job;
                temp_job_machine.machine = i;
                g_array_append_val(rcl, temp_job_machine);
            }

            //}
        }

        /** Choose uniformaly an assignment of a job to a machine */
        int a = g_rand_int_range(rand_, 0, rcl->len);
        int job = g_array_index(rcl, pair_job_machine, a).job;
        int machine = g_array_index(rcl, pair_job_machine, a).machine;
        partlist_insert(temp + machine, new_sol->vlist, jobarray + job);
        g_queue_pop_nth(to_do_list, g_queue_index(to_do_list, jobarray + job));
        g_array_free(rcl, TRUE);
    }

    g_queue_free(to_do_list);
    return 0;
}

int random_assignment(Job *jobarray, int njobs, int nmachines,
                      solution *new_sol, GRand *rand_) {
    int i, val = 0;
    double n;
    partlist *temp = (partlist *) NULL;
    Job *j = (Job *) NULL;
    GQueue *queue = (GQueue *) NULL;
    queue = g_queue_new();

    for (i = 0; i < nmachines; ++i) {
        g_queue_push_head(queue, new_sol->part + i);
    }

    for (i = 0; i < njobs; ++i) {
        j = jobarray + i;
        n = g_rand_double_range(rand_, 0.0, 1.0);

        if (n < 0.8) {
            temp = (partlist *) g_queue_pop_head(queue);
        } else if (n >= 0.8 && n < 0.95) {
            temp = (partlist *) g_queue_pop_nth(queue, 1);
        } else {
            temp = (partlist *)g_queue_pop_nth(queue, 2);
        }

        val = partlist_insert(temp, new_sol->vlist, j);
        CCcheck_val_2(val, "Failed in partlist_insert_order");
        g_queue_insert_sorted(queue, temp, compare_func1, NULL);
    }

CLEAN:
    g_queue_free(queue);
    return val;
}


static int solution_set_c(solution *sol) {
    int val = 0;
    partlist *tmp = (partlist *) NULL;
    Job *j = (Job *) NULL;
    BinomialHeap *heap = binomial_heap_new(BINOMIAL_HEAP_TYPE_MIN,
                                           compare_completion_time);
    CCcheck_NULL_2(heap, "Failed to allocate memory to heap");
    sol->tw = 0;
    sol->b = 0;

    for (unsigned i = 0; i < sol->nmachines; ++i) {
        sol->part[i].c = 0;
        sol->part[i].tw = 0;
        sol->part[i].key = i;
        binomial_heap_insert(heap, sol->part + i);
    }

    for (unsigned i = 0; i < sol->njobs; ++i) {
        j = sol->perm[i];
        tmp = (partlist *) binomial_heap_pop(heap);
        j->index = tmp->machine->len;
        g_ptr_array_add(tmp->machine, j);
        tmp->c += j->processingime;
        sol->c[j->job] = tmp->c;
        tmp->tw += j->weight * (CC_MAX(0, tmp->c - j->duetime));
        sol->tw += j->weight * (CC_MAX(0, tmp->c - j->duetime));
        sol->b += j->duetime * (sol->njobs - i);
        binomial_heap_insert(heap, tmp);
    }


CLEAN:

    if (val) {
        solution_free(sol);
    }

    binomial_heap_free(heap);
    return val;
}

int construct_spt(wctproblem *prob, solution *sol) {
    int val = 0;

    for (unsigned i = 0; i < prob->njobs; ++i) {
        sol->perm[i] = prob->ojobarray[i];
    }

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

    for (unsigned i = 0; i < prob->njobs; ++i) {
        sol->perm[i] = prob->ojobarray[i];
    }

    sol->njobs = prob->njobs;
    sol->nmachines = prob->nmachines;

    val = solution_set_c(sol);
    CCcheck_val_2(val, "failed in solution_set_c");

CLEAN:
    return val;
}


void permutation_solution(GRand *rand_uniform, solution *sol) {
    unsigned i;
    Job *tmp = (Job *) NULL;

    for (i = 0; i <= sol->njobs - 2 ; i++) {
        unsigned j = g_rand_int_range(rand_uniform, 0,
                                      sol->njobs - i);
        CC_SWAP(sol->perm[i], sol->perm[i + j],
                tmp);
    }
}

void solution_forward_insertion(solution *sol, int i, int j) {
    Job *tmp = sol->perm[i];
    memcpy(sol->perm + i, sol->perm + i + 1, (j - i)*sizeof(Job *));
    sol->perm[j] = tmp;
}

void solution_forward_insertion_inverse(solution *sol, int i, int j) {
    Job *tmp = sol->perm[j];
    memcpy(sol->perm + i + 1, sol->perm + i, (j - i)*sizeof(Job *));
    sol->perm[i] = tmp;
}

void local_search_gpi(solution *sol, int *iteration) {
    int move = 1;
    int tw2;
    int b2;
    Job *tmp = (Job *) NULL;


    while (move) {
        tw2 = sol->tw;
        b2 = sol->b;
        move = 0;

        for (unsigned i = 0; i < sol->njobs - 1 && !move; ++i) {
            for (unsigned j = i + 1; j < sol->njobs && !move ; ++j) {
                // if (!move && j - i > 1) {
                //     solution_forward_insertion_inverse(sol, i, j);
                //     solution_set_c(sol);

                //     if (sol->tw < tw2 || (sol->tw == tw2 && sol->b < b2)) {
                //         tw2 = sol->tw;
                //         b2 = sol->b;
                //         move = 1;
                //     } else {
                //         solution_forward_insertion(sol, i, j);
                //     }
                // }
                
                if (!move) {
                    CC_SWAP(sol->perm[i], sol->perm[j], tmp);
                    solution_set_c(sol);

                    if ((sol->tw < tw2) || (sol->tw == tw2 && sol->b < b2)) {
                        tw2 = sol->tw;
                        b2 = sol->b;
                        move = 1;
                    } else {
                        CC_SWAP(sol->perm[i], sol->perm[j], tmp);
                    }
                }

                // if (!move && j - i > 1) {
                //     solution_forward_insertion(sol, i, j);
                //     solution_set_c(sol);

                //     if (sol->tw < tw2 || (sol->tw == tw2 && sol->b < b2)) {
                //         tw2 = sol->tw;
                //         b2 = sol->b;
                //         move = 1;
                //     } else {
                //         solution_forward_insertion_inverse(sol, i, j);
                //     }
                // }





            }
        }

        (*iteration)++;
    }

    printf("number of iterations %d\n", *iteration);
}

int heuristic_rpup(wctproblem *prob) {
    int  val = 0;
    GRand *rand_uniform = g_rand_new_with_seed(1984);
    int N = 4* prob->njobs;
    int i = 1;
    solution *opt_sol = (solution *) NULL;
    solution *sol = (solution *) NULL;
    CCutil_timer test;
    CCutil_init_timer(&test, (char *) NULL);
    int k, l;
    Job *tmp;
    local_search_data *data;

    sol = solution_alloc(prob->nmachines, prob->njobs);
    CCcheck_NULL_2(sol, "Failed to allocate memory");
    val = construct_edd(prob, sol);
    CCcheck_val_2(val, "Failed construct edd");
    data = local_search_data_init(sol);
    for(unsigned i = 0; i < 10; ++i) {
        local_search_swap_inter(sol, data, 1,2);
    }
    local_search_data_free(data);


    prob->opt_sol = solution_alloc(prob->nmachines, prob->njobs);
    CCcheck_NULL_2(prob->opt_sol, "Failed to allocate memory");
    opt_sol = prob->opt_sol;


    CCcheck_val_2(val, "Failed in construct_edd");
    solution_update(opt_sol, sol);

    // while (i < N) {
    //     if (i % 5 == 0) {
    //         permutation_solution(rand_uniform, sol);
    //         solution_set_c(sol);
    //     }

    //     CCutil_start_timer(&test);
    //     local_search_gpi(sol, &i);
    //     CCutil_suspend_timer(&test);
    //     printf("%f\n", test.cum_zeit);

    //     if (sol->tw < opt_sol->tw) {
    //         solution_update(opt_sol, sol);
    //     }

    //     if (i % 20 == 0) {
    //         printf("test iteration %d\n", i);
    //     }

    //     for (unsigned j = 0; j < 3; ++j) {
    //         k = g_rand_int_range(rand_uniform, 0, sol->njobs);

    //         do l = g_rand_int_range(rand_uniform, 0, sol->njobs); while (k == l);

    //         CC_SWAP(sol->perm[k], sol->perm[l], tmp);
    //     }

    //     solution_set_c(sol);

    //     i++;
    // }

    solution_print(prob->opt_sol);

CLEAN:
    solution_free(sol);
    CC_IFFREE(sol, solution);
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
                                      problem->global_lower_bound) / (problem->global_lower_bound);
        partlist_to_Scheduleset(new_sol->part, new_sol->nmachines, new_sol->njobs,
                                &(problem->bestschedule), &(problem->nbestschedule));
    }

    if (problem->global_upper_bound == problem->global_lower_bound) {
        problem->status = optimal;
    }
}

static int add_feasible_solution(wctproblem *problem, solution *new_sol) {
    int val = 0;
    wctdata *root_pd = &(problem->root_pd);
    solution_calc(new_sol, root_pd->jobarray);
    solution_unique(new_sol);


    update_bestschedule(problem, new_sol);

    if (root_pd->ccount == 0 && problem->parms.construct != 0) {
        update_Schedulesets(&root_pd->cclasses, &root_pd->ccount, problem->bestschedule,
                            problem->nbestschedule);
        root_pd->gallocated = root_pd->ccount;
    } else if (problem->parms.construct != 0) {
        partlist_to_Scheduleset(new_sol->part, new_sol->nmachines, new_sol->njobs,
                                &(root_pd->newsets), &(root_pd->nnewsets));
        add_newsets(root_pd);
    }


    return val;
}


int construct_feasible_solutions(wctproblem *problem) {
    int val = 0;
    int iterations = 0;
    wctdata *pd = &(problem->root_pd);
    wctparms *parms = &(problem->parms) ;
    CCutil_timer *timer = &(problem->tot_scatter_search);
    GRand *rand1 = g_rand_new_with_seed(1984);
    GRand *rand2 = g_rand_new_with_seed(1654651);
    CCutil_start_timer(timer);

    //while (1) {
        iterations++;
        solution *new_sol = (solution *) NULL;
        new_sol = solution_alloc(pd->nmachines, pd->njobs);
        CCcheck_NULL(new_sol, "Failed to allocate")

        if (problem->status == no_sol) {

        } else {
            if (g_rand_boolean(rand1)) {
                random_assignment(pd->jobarray, pd->njobs, pd->nmachines, new_sol, rand2);
            } else {
                random_rcl_assignment(pd->jobarray, pd->njobs, pd->nmachines, new_sol, rand2);
            }
        }

        val = add_feasible_solution(problem, new_sol);
        CCcheck_val(val, "Failed in add_feasible_solution");
        //break;
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
