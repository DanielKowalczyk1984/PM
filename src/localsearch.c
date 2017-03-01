#include <assert.h>
#include "wct.h"

int move(Job *j, partlist *m_j, partlist *m_i);

int sort_jobs(gconstpointer a, gconstpointer b) {
    int aa = (((const Job *)a)->job);
    int bb = (((const Job *)b)->job);
    return aa - bb;
}

int move(Job *j, partlist *m_j, partlist *m_i) {
    int nb_job = j->job;
    return 0;
}

int k_l_move_general(Job **K_jobs, Job **L_jobs, partlist *m_k, partlist *m_l,
                     solution *sol, int k, int l) {
    int i, val = 0;
    Job *ptr_job = (Job *) NULL;
    Job *ptr_job2 = (Job *) NULL;
    GList *list = (GList *) NULL;
    GList *it = (GList *) NULL;
    GList *it2 = (GList *) NULL;

    if (l > 0) {
        m_l = sol->vlist[L_jobs[0]->job].part;
    }

    for (i = 0;  i < k; i++) {
        list = g_list_insert_sorted(list, K_jobs[i], sort_jobs);
    }

    if (l > 0) {
        for (i = 0; i < l; i++) {
            list = g_list_insert_sorted(list, L_jobs[i], sort_jobs);
        }
    }

    for (it = list; it; it = it->next) {
        ptr_job = (Job *) it->data;

        if (l > 0) {
            if (sol->vlist[ptr_job->job].part == m_k) {
                val += move(ptr_job, m_k, m_l);
            } else {
                val += move(ptr_job, m_l, m_k);
            }

            it2 = list;

            while (it2 != it) {
                ptr_job2 = (Job *)it2->data;

                if (sol->vlist[ptr_job->job].part != sol->vlist[ptr_job2->job].part) {
                    val += 2 * ptr_job->weight * ptr_job2->processingime;
                } else {
                    val -= 2 * ptr_job->weight * ptr_job2->processingime;
                }

                it2 = it2->next;
            }
        } else {
            val += move(ptr_job, m_k, m_l);
            it2 = list;

            while (it2 != it) {
                ptr_job2 = (Job *)it2->data;

                if (sol->vlist[ptr_job->job].part != sol->vlist[ptr_job2->job].part) {
                    val += 2 * ptr_job->weight * ptr_job2->processingime;
                } else {
                    val -= 2 * ptr_job->weight * ptr_job2->processingime;
                }

                it2 = it2->next;
            }
        }
    }

    g_list_free(list);
    return val;
}

int local_search_machine_general_best(solution *sol, int lowerbound, int k,
                                      int l) {
    int i, j, n1, n2, it1, it2, K_flag, L_flag, moved, val = 0;
    int nbiter = 0;
    int njobs = sol->njobs;
    int nmachines = sol->nmachines;
    partlist *L_max = (partlist *) NULL;
    partlist *K_max = (partlist *) NULL;
    CCutil_timer time_move;
    CCutil_init_timer(&time_move, NULL);
    int max;
    int improvement = 0;
    Job **K_jobs = (Job **) NULL;
    Job **L_jobs = (Job **) NULL;
    Job **K_jobs_max = (Job **) NULL;
    Job **L_jobs_max = (Job **) NULL;
    int *K_subset = (int *) NULL;
    int *L_subset = (int *) NULL;
    partlist *temp_m1, *temp_m2;
    K_subset = CC_SAFE_MALLOC(k + 1, int);
    K_jobs = CC_SAFE_MALLOC(k, Job *);
    K_jobs_max = CC_SAFE_MALLOC(k, Job *);

    if (l > 0) {
        L_jobs = CC_SAFE_MALLOC(l, Job *);
        L_subset = CC_SAFE_MALLOC(l + 1, int);
        L_jobs_max = CC_SAFE_MALLOC(l, Job *);
    }

    if (dbg_lvl() > 0) {
        printf("Executing local search with k = %d and l = %d\n", k, l);
    }

    moved = 1;
    CCutil_start_timer(&time_move);

    while (moved && sol->tw != lowerbound) {
        nbiter++;
        moved = 0;
        max = 0;
        int temp;

        for (i = 0; i < nmachines - 1; ++i) {
            for (j = i + 1; j < nmachines; ++j) {
                temp_m1 = sol->part + i;
                temp_m2 = sol->part + j;
                n1 = g_queue_get_length(temp_m1->list);

                if (k > n1) {
                    continue;
                }

                k_subset_init(n1, k, K_subset, &K_flag);

                while (K_flag) {
                    for (it1 = 0; it1 < k; ++it1) {
                        K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                    }

                    if (l == 0) {
                        temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                        if (temp > max) {
                            for (it1 = 0; it1 < k; ++it1) {
                                K_jobs_max[it1] = K_jobs[it1];
                            }

                            L_max = temp_m2;
                            moved = 1;
                            max = temp;
                        }
                    } else if (l > 0) {
                        n2 = g_queue_get_length(temp_m2->list);

                        if (l > n2) {
                            continue;
                        }

                        k_subset_init(n2, l, L_subset, &L_flag);

                        while (L_flag) {
                            for (it2 = 0; it2 < l; ++it2) {
                                L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                            }

                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                for (it1 = 0; it1 < l; ++it1) {
                                    L_jobs_max[it1] = L_jobs[it1];
                                }

                                L_max = temp_m2;
                                K_max = temp_m1;
                                moved = 1;
                                max = temp;
                            }

                            k_subset_lex_successor(n2, l, L_subset, &L_flag);
                        }
                    }

                    k_subset_lex_successor(n1, k, K_subset, &K_flag);
                }

                temp_m1 = sol->part + j;
                temp_m2 = sol->part + i;
                n1 = g_queue_get_length(temp_m1->list);

                if (k > n1) {
                    continue;
                }

                k_subset_init(n1, k, K_subset, &K_flag);

                while (K_flag) {
                    for (it1 = 0; it1 < k; ++it1) {
                        K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                    }

                    if (l == 0) {
                        temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                        if (temp > max) {
                            for (it1 = 0; it1 < k; ++it1) {
                                K_jobs_max[it1] = K_jobs[it1];
                            }

                            L_max = temp_m2;
                            moved = 1;
                            max = temp;
                        }
                    } else if (l > 0) {
                        n2 = g_queue_get_length(temp_m2->list);

                        if (l > n2) {
                            continue;
                        }

                        k_subset_init(n2, l, L_subset, &L_flag);

                        while (L_flag) {
                            for (it2 = 0; it2 < l; ++it2) {
                                L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                            }

                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                for (it1 = 0; it1 < l; ++it1) {
                                    L_jobs_max[it1] = L_jobs[it1];
                                }

                                L_max = temp_m2;
                                K_max = temp_m1;
                                max = temp;
                                moved = 1;
                            }

                            k_subset_lex_successor(n2, l, L_subset, &L_flag);
                        }
                    }

                    k_subset_lex_successor(n1, k, K_subset, &K_flag);
                }
            }
        }

        if (moved && max > 0) {
            improvement += max;

            if (l == 0) {
                for (it1 = 0; it1 < k; ++it1) {
                    //partlist_move_order(L_max, sol->vlist, K_jobs_max[it1], njobs);
                }
            } else if (l > 0) {
                for (it1 = 0; it1 < k; ++it1) {
                    //partlist_move_order(L_max, sol->vlist, K_jobs_max[it1], njobs);
                }

                for (it1 = 0; it1 < l; ++it1) {
                    //partlist_move_order(K_max, sol->vlist, L_jobs_max[it1], njobs);
                }
            }

            sol->tw -= max;
        }
    };

    CCutil_stop_timer(&time_move, 0);

    if (dbg_lvl()) {
        printf("local search %d  - %d -> number of iterations = %d, objective = %d, time = %f, time per iteration = %f, improvement = %d\n",
               k, l, nbiter, sol->tw, time_move.cum_zeit,
               time_move.cum_zeit / nbiter, improvement);
    }

    CC_IFFREE(L_jobs, Job *);
    CC_IFFREE(K_jobs, Job *);
    CC_IFFREE(K_jobs_max, Job *);
    CC_IFFREE(L_jobs_max, Job *);
    CC_IFFREE(L_subset, int);
    CC_IFFREE(K_subset, int);
    return val;
}

int local_search_machine_general_first(solution *sol, int lowerbound, int k,
                                       int l) {
    int i, j, n1, n2, it1, it2, K_flag, L_flag, moved, val = 0;
    int nbiter = 0;
    int njobs = sol->njobs;
    int nmachines = sol->nmachines;
    partlist *L_max = (partlist *) NULL;
    partlist *K_max = (partlist *) NULL;
    int max;
    CCutil_timer time_move;
    CCutil_init_timer(&time_move, NULL);
    int improvement = 0;
    Job **K_jobs = (Job **) NULL;
    Job **L_jobs = (Job **) NULL;
    Job **K_jobs_max = (Job **) NULL;
    Job **L_jobs_max = (Job **) NULL;
    int *K_subset = (int *) NULL;
    int *L_subset = (int *) NULL;
    partlist *temp_m1, *temp_m2;
    K_subset = CC_SAFE_MALLOC(k + 1, int);
    K_jobs = CC_SAFE_MALLOC(k, Job *);
    K_jobs_max = CC_SAFE_MALLOC(k, Job *);

    if (l > 0) {
        L_jobs = CC_SAFE_MALLOC(l, Job *);
        L_subset = CC_SAFE_MALLOC(l + 1, int);
        L_jobs_max = CC_SAFE_MALLOC(l, Job *);
    }

    if (dbg_lvl()) {
        printf("Executing local search with k = %d and l = %d\n", k, l);
    }

    moved = 1;
    CCutil_start_timer(&time_move);

    while (moved && sol->tw != lowerbound) {
        nbiter++;
        moved = 0;
        max = 0;
        int temp;

        for (i = 0; i < nmachines - 1 && !moved; ++i) {
            for (j = i + 1; j < nmachines && !moved; ++j) {
                temp_m1 = sol->part + i;
                temp_m2 = sol->part + j;
                n1 = g_queue_get_length(temp_m1->list);

                if (k > n1) {
                    continue;
                }

                k_subset_init(n1, k, K_subset, &K_flag);

                while (K_flag && !moved) {
                    for (it1 = 0; it1 < k; ++it1) {
                        K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                    }

                    if (l == 0) {
                        temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                        if (temp > max) {
                            for (it1 = 0; it1 < k; ++it1) {
                                K_jobs_max[it1] = K_jobs[it1];
                            }

                            L_max = temp_m2;
                            moved = 1;
                            max = temp;
                        }
                    } else if (l > 0) {
                        n2 = g_queue_get_length(temp_m2->list);
                        k_subset_init(n2, l, L_subset, &L_flag);

                        if (l > n2) {
                            K_flag = 0;
                            continue;
                        };

                        while (L_flag && !moved) {
                            for (it2 = 0; it2 < l; ++it2) {
                                L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                            }

                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                for (it1 = 0; it1 < l; ++it1) {
                                    L_jobs_max[it1] = L_jobs[it1];
                                }

                                L_max = temp_m2;
                                K_max = temp_m1;
                                moved = 1;
                                max = temp;
                            }

                            k_subset_lex_successor(n2, l, L_subset, &L_flag);
                        }
                    }

                    k_subset_lex_successor(n1, k, K_subset, &K_flag);
                }

                if (k != l) {
                    temp_m1 = sol->part + j;
                    temp_m2 = sol->part + i;
                    n1 = g_queue_get_length(temp_m1->list);

                    if (k > n1) {
                        continue;
                    }

                    k_subset_init(n1, k, K_subset, &K_flag);

                    while (K_flag && !moved) {
                        for (it1 = 0; it1 < k; ++it1) {
                            K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                        }

                        if (l == 0) {
                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                L_max = temp_m2;
                                moved = 1;
                                max = temp;
                            }
                        } else if (l > 0) {
                            n2 = g_queue_get_length(temp_m2->list);

                            if (l > n2) {
                                K_flag = 0;
                                continue;
                            };

                            k_subset_init(n2, l, L_subset, &L_flag);

                            while (L_flag && !moved) {
                                for (it2 = 0; it2 < l; ++it2) {
                                    L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                                }

                                temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                                if (temp > max) {
                                    for (it1 = 0; it1 < k; ++it1) {
                                        K_jobs_max[it1] = K_jobs[it1];
                                    }

                                    for (it1 = 0; it1 < l; ++it1) {
                                        L_jobs_max[it1] = L_jobs[it1];
                                    }

                                    L_max = temp_m2;
                                    K_max = temp_m1;
                                    max = temp;
                                    moved = 1;
                                }

                                k_subset_lex_successor(n2, l, L_subset, &L_flag);
                            }
                        }

                        k_subset_lex_successor(n1, k, K_subset, &K_flag);
                    }
                }
            }
        }

        if (moved && max > 0) {
            improvement += max;

            for (it1 = 0; it1 < k; ++it1) {
                //partlist_move_order(L_max, sol->vlist, K_jobs_max[it1], njobs);
            }

            if (l > 0) {
                for (it1 = 0; it1 < l; ++it1) {
                    //partlist_move_order(K_max, sol->vlist, L_jobs_max[it1], njobs);
                }
            }

            sol->tw -= max;
        }
    };

    CCutil_stop_timer(&time_move, 0);

    if (dbg_lvl()) {
        printf("local search %d  - %d -> number of iterations = %d, objective = %d, time = %f, time per iteration = %f, improvement = %d\n",
               k, l, nbiter, sol->tw, time_move.cum_zeit,
               time_move.cum_zeit / nbiter, improvement);
    }

    CC_IFFREE(L_jobs, Job *);
    CC_IFFREE(K_jobs, Job *);
    CC_IFFREE(K_jobs_max, Job *);
    CC_IFFREE(L_jobs_max, Job *);
    CC_IFFREE(L_subset, int);
    CC_IFFREE(K_subset, int);
    return val;
}

void localsearch_wrap(solution *sol, int lowerbound, int best) {
    if (best) {
        local_search_machine_general_best(sol, lowerbound, 1, 0);
        local_search_machine_general_best(sol, lowerbound, 1, 1);
        local_search_machine_general_best(sol, lowerbound, 2, 0);
        local_search_machine_general_best(sol, lowerbound, 2, 1);
    } else {
        local_search_machine_general_first(sol, lowerbound, 1, 0);
        local_search_machine_general_first(sol, lowerbound, 1, 1);
        local_search_machine_general_first(sol, lowerbound, 2, 0);
        local_search_machine_general_first(sol, lowerbound, 2, 1);
    }
}

void localsearch_random_k(solution *sol, int lowerbound, int nb) {
    int i, j, l, k, tot;
    int **matrix = (int **) NULL;
    matrix = CC_SAFE_MALLOC(nb, int *);

    for (i = 0; i < nb; ++i) {
        matrix[i] = CC_SAFE_MALLOC(nb, int);

        for (j = 0; j < nb; ++j) {
            matrix[i][j] = 0;
        }
    }

    tot = ((2 + nb) * (nb - 1)) / 2 - 1;
    k = g_random_int_range(1, nb);
    l = g_random_int_range(0, k + 1);
    local_search_machine_general_first(sol, lowerbound, k, l);
    matrix[k][l] = 1;

    for (i = 0; i < tot && !(sol->tw == lowerbound); ++i) {
        while (matrix[k][l] != 0) {
            k = g_random_int_range(1, nb);
            l = g_random_int_range(0, k + 1);
        }

        local_search_machine_general_first(sol, lowerbound, k, l);
        matrix[k][l] = 1;
    }

    for (i = 0; i < nb; ++i) {
        CC_IFFREE(matrix[i], int);
    }

    CC_IFFREE(matrix, int *);
}

local_search_data *local_search_data_init(solution *sol) {
    int val = 0;
    local_search_data *data;
    int nmachines = sol->nmachines;
    int i, j;

    data = CC_SAFE_MALLOC(1, local_search_data);
    CCcheck_NULL_2(data, "Failed to allocate memory");
    data->nmachines = nmachines;

    data->W = CC_SAFE_MALLOC(nmachines, int *);
    data->g = CC_SAFE_MALLOC(nmachines, GList **);
    data->processing_list = CC_SAFE_MALLOC(nmachines, processing_list_data *);
    data->len_list_slope_t = 0;
    data->list_slope_t = CC_SAFE_MALLOC((sol->njobs + 1) * (sol->njobs + 1),
                                        slope_t);
    data->njobs = CC_SAFE_MALLOC(nmachines, int);

    for (i = 0; i < nmachines; ++i) {
        data->njobs[i] = sol->part[i].machine->len;
        data->W[i] = CC_SAFE_MALLOC(data->njobs[i], int);
        data->g[i] = CC_SAFE_MALLOC(data->njobs[i], GList *);
        data->processing_list[i] = (processing_list_data *) NULL;

        for (j = 0; j < sol->part[i].machine->len; ++j) {
            data->g[i][j] = (GList *) NULL;
        }
    }


CLEAN:

    if (val) {
        for (i = 0; i < nmachines; ++i) {
            for (j = 0; j < sol->part[i].machine->len; ++j) {
                g_list_free(data->g[i][j]);
            }

            CC_IFFREE(data->g[i], GList *);
            CC_IFFREE(data->W[i], int);
            if(data->processing_list[i] != NULL) {
                CC_IFFREE(data->processing_list[i], processing_list_data);
            }
        }

        CC_IFFREE(data->g, GList **);
        CC_IFFREE(data->W, int *);
        CC_IFFREE(data->processing_list, processing_list_data *);
        CC_IFFREE(data->njobs, int);
        CC_IFFREE(data->list_slope_t, slope_t);
        CC_IFFREE(data, local_search_data);
    }

    return data;
}

void local_search_data_free(local_search_data *data) {
    int nmachines, i, j;

    if (data != (local_search_data *) NULL) {
        nmachines = data->nmachines;

        for (i = 0; i < nmachines; ++i) {
            for (j = 0; j < data->njobs[i]; ++j) {
                g_list_free(data->g[i][j]);
            }

            CC_IFFREE(data->g[i], GList *);
            CC_IFFREE(data->W[i], int);
            if(data->processing_list[i] != NULL) {
                CC_IFFREE(data->processing_list[i], processing_list_data);
            }
        }

        CC_IFFREE(data->g, GList **);
        CC_IFFREE(data->W, int *);
        CC_IFFREE(data->processing_list, processing_list_data *);
        CC_IFFREE(data->njobs, int);
        CC_IFFREE(data->list_slope_t, slope_t);
        CC_IFFREE(data, local_search_data);
    }
}



int local_search_data_calculate(solution *sol, local_search_data *data) {
    int val = 0;
    int nmachines;
    Job *tmp;

    if (sol == NULL || data == NULL || sol->nmachines != data->nmachines) {
        val = 1;
        printf("Not compatible data structures\n");
        return val;
    }

    nmachines = data->nmachines;

    /**
     * Create W for every machine and job
     */
    for (unsigned i = 0; i < nmachines; ++i) {
        tmp = (Job *)g_ptr_array_index(sol->part[i].machine, 0);
        data->W[i][0] = tmp->weight * CC_MAX(sol->c[tmp->job] - tmp->duetime, 0);

        for (unsigned j = 1; j < data->njobs[i]; j++) {
            tmp = (Job *)g_ptr_array_index(sol->part[i].machine, j);
            data->W[i][j] = data->W[i][j - 1] + tmp->weight * CC_MAX(
                                sol->c[tmp->job] - tmp->duetime, 0);
        }
    }


    return val;
}

int compare_process_list(gconstpointer a, gconstpointer b) {
    processing_list_data *x = (processing_list_data *) a;
    processing_list_data *y = (processing_list_data *) b;

    if (x->p > y->p) {
        return -1;
    } else if (x->p < y->p) {
        return 1;
    }

    return 0;

}

int local_search_create_processing_list(solution *sol, local_search_data *data,
                                        int l) {
    int val = 0;
    Job *j1, *j2;
    int C;

    for (unsigned i = 0; i < data->nmachines; ++i) {
        C = 0;

        for (unsigned j = 0; j < l; ++j) {
            j1 = (Job *) g_ptr_array_index(sol->part[i].machine, j);
            C += j1->processingime;
        }

        data->processing_list[i] = CC_SAFE_MALLOC(data->njobs[i] - l, processing_list_data);

        for (unsigned j = 0; j < data->njobs[i] - l ; ++j) {
            j1 = (Job *) g_ptr_array_index(sol->part[i].machine, j);
            j2 = (Job *) g_ptr_array_index(sol->part[i].machine, j + l);
            data->processing_list[i][j].pos = j;
            data->processing_list[i][j].p = C;
            C = C -  j1->processingime + j2->processingime;
        }

        qsort(data->processing_list[i], data->njobs[i] - l, sizeof(processing_list_data), compare_process_list);
        for(unsigned j = 0; j < data->njobs[i]-l; ++j) {
            printf("%d  \n",data->processing_list[i][j].p);
        }
    }

    return val;
}

int local_search_compare_lateness(gconstpointer a, gconstpointer b,
                                  gpointer data) {
    int *data_x = (int *) data;
    Job *x = *(Job **) a;
    Job *y = *(Job **) b;

    if (data_x[x->job] - x->duetime < data_x[y->job] - y->duetime) {
        return 1;
    } else if (data_x[x->job] - x->duetime > data_x[y->job] - y->duetime) {
        return -1;
    }

    return 0;
}

static void local_search_add_slope_t(local_search_data *data, int b1, int b2, int c, int alpha, int i, int j){
    data->list_slope_t[data->len_list_slope_t].alpha = alpha;
    data->list_slope_t[data->len_list_slope_t].c = c;
    data->list_slope_t[data->len_list_slope_t].b1 = b1;
    data->list_slope_t[data->len_list_slope_t].b2 = b2;
    data->g[i][j] = g_list_append(data->g[i][j], data->list_slope_t + data->len_list_slope_t++);
}

int local_search_create_g(solution *sol, local_search_data *data) {
    int val = 0;
    int nmachines = sol->nmachines;
    Job *tmp;
    int t1, t2;
    int tw;
    int w;

    for (unsigned i = 0; i < nmachines; ++i) {
        int n_k = sol->part[i].machine->len;
        GPtrArray *lateness_sort =  g_ptr_array_new();

        for (unsigned j = 0; j < n_k; ++j) {
            g_ptr_array_add(lateness_sort, g_ptr_array_index(sol->part[i].machine, j));
        }

        g_ptr_array_sort_with_data(lateness_sort, local_search_compare_lateness,
                                   sol->c);

        int P = 0;

        for (unsigned j = 0; j < n_k ; ++j) {
            tw = 0;
            w = 0;
            t1 = 0;
            int k;
            int move;

            move = 1;
            tmp = (Job *) g_ptr_array_index(lateness_sort, j);
            for(k = j  ; k < n_k && move; ) {
                move = tmp->weight * (sol->c[tmp->job] - P - tmp->duetime) > 0;
                if(move) {
                    tw += tmp->weight * (sol->c[tmp->job] - P - tmp->duetime);
                    w += tmp->weight;
                    k++;
                    tmp = (Job *) g_ptr_array_index(lateness_sort, k);
                } 
            }


            t2 = tmp->duetime - sol->c[tmp->job] + P;
            local_search_add_slope_t(data, t1, t2, tw, w, i, j);


            for (unsigned l = k ; l < n_k; ) {
                tw = tw + w * (t2 - t1);
                t1 = t2;

                move = 1;
                tmp = (Job *) g_ptr_array_index(lateness_sort, l);
                while(move) {
                    w += tmp->weight;
                    l++;
                    if(l == n_k) {
                        move = 0;
                        t2 = INT_MAX;
                    } else {
                        tmp = (Job *) g_ptr_array_index(lateness_sort, l);
                        t2 = tmp->duetime - sol->c[tmp->job] + P;
                        move = (t1 == t2); 
                    }
                }
  
                local_search_add_slope_t(data, t1, t2, tw, w, i, j);
            }

            tmp = (Job *) g_ptr_array_index(sol->part[i].machine, j);
            P += tmp->processingime;

        }
        g_ptr_array_free(lateness_sort, TRUE);
    }

    return val;
}

static int compute_g(slope_t *x, int t){
    return x->c + x->alpha*(t - x->b1);
}

void local_search_forward_insertion(solution *sol, local_search_data *data, int l){
    int pos,p;
    slope_t *tmp;
    int **g,**h, *gg, **hh;
    int move;
    GList **iterators;
    Job *tmp_j;

    g = CC_SAFE_MALLOC(sol->njobs, int*);
    h = CC_SAFE_MALLOC(sol->njobs, int*);
    hh = CC_SAFE_MALLOC(sol->njobs, int*);
    iterators = CC_SAFE_MALLOC(sol->njobs, GList*);
    gg = CC_SAFE_MALLOC(sol->njobs, int);
    for(unsigned i = 0; i < sol->njobs; ++i) {
        g[i] = CC_SAFE_MALLOC(sol->njobs, int);
        h[i] = CC_SAFE_MALLOC(sol->njobs, int);
        hh[i] = CC_SAFE_MALLOC(sol->njobs, int);
    }

 

    //for(unsigned k = 0; k < sol->nmachines; ++k) {
        /** compute g */
        for(unsigned i = 0; i < data->njobs[0] - l; ++i) {
            pos = data->processing_list[0][i].pos;
            p = data->processing_list[0][i].p;
            GList *it = data->g[0][pos];
            for(unsigned j = pos + l; j < data->njobs[0]; ++j) {
                tmp = (slope_t*) it->data;
                tmp_j = (Job *) g_ptr_array_index(sol->part[0].machine, j);
                move = !(tmp->b1 <= sol->c[tmp_j->job] - p && tmp->b2 >= sol->c[tmp_j->job] - p);
                while(move){
                    it = it->next;
                    if(it != (GList *) NULL) {
                        tmp = (slope_t *) it->data;
                        move = !(tmp->b1 <= sol->c[tmp_j->job] - p && tmp->b2 >= sol->c[tmp_j->job] - p);
                    } else {
                        move = 0;
                    }
                }
                g[pos][j] = compute_g(tmp,sol->c[tmp_j->job] - p);

            }
        }
        /** compute h */
        for(unsigned i = 0; i < data->njobs[0] - l; ++i) {
            GList *it = data->g[0][i + l];
            for(unsigned j = i + l; j < data->njobs[0]; ++j) {
                tmp = (slope_t*) it->data;
                tmp_j = (Job *) g_ptr_array_index(sol->part[0].machine, j);
                move = !(tmp->b1 <= sol->c[tmp_j->job] && tmp->b2 >= sol->c[tmp_j->job]);
                while(move){
                    it = it->next;
                    if(it != (GList *) NULL) {
                        tmp = (slope_t *) it->data;
                        move = !(tmp->b1 <= sol->c[tmp_j->job] && tmp->b2 >= sol->c[tmp_j->job]);
                    } else {
                        move = 0;
                    }
                }
                h[i][j] = compute_g(tmp,sol->c[tmp_j->job]);

            }
        }

        /** compute gg */
        for(unsigned i = 0; i < data->njobs[0] - l; ++i) {
            GList *it = data->g[0][i + l];
            tmp = (slope_t*) it->data;
            int c;
            if(i != 0) {
                tmp_j = (Job *) g_ptr_array_index(sol->part[0].machine, i - 1);
                c  = sol->c[tmp_j->job];
            } else {
                c = 0;
            }
            move = !(tmp->b1 <= c && tmp->b2 >= c);
            while(move){
                it = it->next;
                if(it != (GList *) NULL) {
                    tmp = (slope_t *) it->data;
                    move = !(tmp->b1 <= c && tmp->b2 >= c);
                } else {
                    move = 0;
                }
            }
            gg[i] = compute_g(tmp,c);

        }

        /** compute hh */
        for(unsigned j = l; j < data->njobs[0] - 1; ++j) {
            iterators[j] = data->g[0][j + 1];
        }
        iterators[data->njobs[0] - 1] = (GList *) NULL;

        for(unsigned i = 0; i < data->njobs[0] - l; ++i) {
            pos = data->processing_list[0][i].pos;
            p = data->processing_list[0][i].p;
            for(unsigned j = pos + l; j < data->njobs[0]; ++j) {
                if(iterators[j] != (GList *) NULL) {
                    tmp = (slope_t*) iterators[j]->data;
                    tmp_j = (Job *) g_ptr_array_index(sol->part[0].machine, j);
                    move = !(tmp->b1 <= sol->c[tmp_j->job] - p && tmp->b2 >= sol->c[tmp_j->job] - p);
                    while(move){
                        iterators[j] = iterators[j]->next;
                        if(iterators[j] != (GList *) NULL) {
                            tmp = (slope_t *) iterators[j]->data;
                            move = !(tmp->b1 <= sol->c[tmp_j->job] - p && tmp->b2 >= sol->c[tmp_j->job] - p);
                        } else {
                            move = 0;
                        }
                    }
                    hh[pos][j] = compute_g(tmp,sol->c[tmp_j->job] - p);
                } else {
                    hh[pos][j] = 0;
                }

            }
        }

    //}

    int max = sol->part[0].tw;
    int i_best,j_best;
    for(unsigned i = 0; i < data->njobs[0] - l; ++i) {
        for(unsigned j = i + l; j < data->njobs[0]; ++j) {
            int tmp = 0;
            if(i != 0) {
                tmp += data->W[0][i - 1];
            }

            tmp += g[i][j] - h[i][j];
            tmp += gg[i] - hh[i][j];
            tmp += data->W[0][data->njobs[0] - 1]  - data->W[0][j- 1];

            printf("test %d %d\n",tmp, sol->part[0].tw);
            if(tmp < max) {
                max = tmp;
                i_best = i;
                j_best = j;
            }
        }
    }

    solution_print(sol);
    printf("best: %d %d %d\n", i_best,j_best, max);
    getchar();

    for(unsigned i = 0; i < sol->njobs; ++i) {
        CC_IFFREE(g[i], int);
        CC_IFFREE(h[i], int);
        CC_IFFREE(hh[i], int);

    }
    CC_IFFREE(g, int*);
    CC_IFFREE(h, int*);
    CC_IFFREE(gg, int);
    CC_IFFREE(hh, int*);
    CC_IFFREE(iterators, GList*);
}
