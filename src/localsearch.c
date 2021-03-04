#include <localsearch.h>
#include <util.h>

int **  g, **h, *gg, **hh;
GList** iterators;

int **B2_1, **B2_2, **B3_1, **B3_2, **B4_1, **B4_2, **B5_1, **B5_2, **B6_1;
int*  B3_1_;

void alloc_all(Solution* sol) {
    g = CC_SAFE_MALLOC(sol->nb_jobs, int*);
    h = CC_SAFE_MALLOC(sol->nb_jobs, int*);
    hh = CC_SAFE_MALLOC(sol->nb_jobs, int*);
    gg = CC_SAFE_MALLOC(sol->nb_jobs, int);
    iterators = CC_SAFE_MALLOC(sol->nb_jobs, GList*);

    for (int i = 0; i < sol->nb_jobs; ++i) {
        g[i] = CC_SAFE_MALLOC(sol->nb_jobs, int);
        h[i] = CC_SAFE_MALLOC(sol->nb_jobs, int);
        hh[i] = CC_SAFE_MALLOC(sol->nb_jobs, int);
    }
    B2_1 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B2_2 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B3_1 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B3_2 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B4_1 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B4_2 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B5_1 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B5_2 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B6_1 = CC_SAFE_MALLOC(sol->nb_jobs + 1, int*);
    B3_1_ = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);

    for (int i = 0; i < sol->nb_jobs + 1; ++i) {
        B2_1[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B2_2[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B3_1[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B3_2[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B4_1[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B4_2[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B5_1[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B5_2[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
        B6_1[i] = CC_SAFE_MALLOC(sol->nb_jobs + 1, int);
    }
}

void free_all(Solution* sol) {
    for (int i = 0; i < sol->nb_jobs; ++i) {
        CC_IFFREE(g[i], int);
        CC_IFFREE(h[i], int);
        CC_IFFREE(hh[i], int);
    }

    CC_IFFREE(g, int*);
    CC_IFFREE(h, int*);
    CC_IFFREE(gg, int);
    CC_IFFREE(hh, int*);
    CC_IFFREE(iterators, GList*);

    for (int i = 0; i < sol->nb_jobs + 1; ++i) {
        CC_IFFREE(B2_1[i], int);
        CC_IFFREE(B2_2[i], int);
        CC_IFFREE(B3_1[i], int);
        CC_IFFREE(B3_2[i], int);
        CC_IFFREE(B4_1[i], int);
        CC_IFFREE(B4_2[i], int);
        CC_IFFREE(B5_1[i], int);
        CC_IFFREE(B5_2[i], int);
        CC_IFFREE(B6_1[i], int);
    }

    CC_IFFREE(B2_1, int*);
    CC_IFFREE(B2_2, int*);
    CC_IFFREE(B3_1, int*);
    CC_IFFREE(B3_2, int*);
    CC_IFFREE(B4_1, int*);
    CC_IFFREE(B4_2, int*);
    CC_IFFREE(B5_1, int*);
    CC_IFFREE(B5_2, int*);
    CC_IFFREE(B6_1, int*);

    CC_IFFREE(B3_1_, int);
}

int compare_process_list(gconstpointer a, gconstpointer b);
int compare_process_list_b(gconstpointer a, gconstpointer b);
int local_search_compare_lateness(gconstpointer a,
                                  gconstpointer b,
                                  gpointer      data);

static void destroy_slope_t(gpointer data) {
    slope_t* tmp = (slope_t*)data;
    CC_IFFREE(tmp, slope_t);
}

static int compute_g(GList** it, int t) {
    slope_t* x = ((*it) != NULL) ? (slope_t*)(*it)->data : NULL;
    return (x != NULL) ? x->c + x->alpha * (t - x->b1) : 0;
}

static void compute_it(GList** it, int c) {
    if ((*it) != NULL) {
        slope_t* tmp = (slope_t*)(*it)->data;
        int      move = !(tmp->b1 <= c && tmp->b2 >= c);

        while (move) {
            *it = (*it)->next;

            if (*it != (GList*)NULL) {
                tmp = (slope_t*)(*it)->data;
                move = !(tmp->b1 <= c && tmp->b2 > c);
            } else {
                move = 0;
            }
        }
    }
}

int compare_process_list(gconstpointer a, gconstpointer b) {
    const processing_list_data* x = (const processing_list_data*)a;
    const processing_list_data* y = (const processing_list_data*)b;

    if (x->p > y->p) {
        return -1;
    } else if (x->p < y->p) {
        return 1;
    }

    return 0;
}

int compare_process_list_b(gconstpointer a, gconstpointer b) {
    const processing_list_data* x = (const processing_list_data*)a;
    const processing_list_data* y = (const processing_list_data*)b;

    if (x->p > y->p) {
        return 1;
    } else if (x->p < y->p) {
        return -1;
    }

    return 0;
}

int local_search_compare_lateness(gconstpointer a,
                                  gconstpointer b,
                                  gpointer      data) {
    int* data_x = (int*)data;
    Job* x = *(Job* const*)a;
    Job* y = *(Job* const*)b;

    if (data_x[x->job] - x->due_time < data_x[y->job] - y->due_time) {
        return 1;
    } else if (data_x[x->job] - x->due_time > data_x[y->job] - y->due_time) {
        return -1;
    }

    return 0;
}

static void local_search_add_slope_t(local_search_data* data,
                                     int                b1,
                                     int                b2,
                                     int                c,
                                     int                alpha,
                                     int                i,
                                     int                j) {
    slope_t* tmp = CC_SAFE_MALLOC(1, slope_t);
    tmp->alpha = alpha;
    tmp->c = c;
    tmp->b1 = b1;
    tmp->b2 = b2;
    data->g[i][j] = g_list_append(data->g[i][j], tmp);
}

local_search_data* local_search_data_init(int njobs, int nb_machines) {
    int                val = 0;
    local_search_data* data = (local_search_data*)NULL;
    int                i, j;
    data = CC_SAFE_MALLOC(1, local_search_data);
    CCcheck_NULL_2(data, "Failed to allocate memory");
    data->nmachines = nb_machines;
    data->W = CC_SAFE_MALLOC(nb_machines, int*);
    data->g = CC_SAFE_MALLOC(nb_machines, GList**);
    data->processing_list_1 =
        CC_SAFE_MALLOC(nb_machines, processing_list_data*);
    data->processing_list_2 =
        CC_SAFE_MALLOC(nb_machines, processing_list_data*);
    data->njobs = njobs;

    for (i = 0; i < nb_machines; ++i) {
        data->W[i] = CC_SAFE_MALLOC(njobs, int);
        data->g[i] = CC_SAFE_MALLOC(njobs, GList*);
        data->processing_list_1[i] =
            CC_SAFE_MALLOC(njobs, processing_list_data);
        data->processing_list_2[i] =
            CC_SAFE_MALLOC(njobs, processing_list_data);

        for (j = 0; j < njobs; ++j) {
            data->g[i][j] = (GList*)NULL;
        }
    }
    data->iterations = 0;

CLEAN:

    if (val && data) {
        for (i = 0; i < nb_machines; ++i) {
            for (j = 0; j < njobs; ++j) {
                g_list_free_full(data->g[i][j], destroy_slope_t);
            }

            CC_IFFREE(data->g[i], GList*);
            CC_IFFREE(data->W[i], int);
            CC_IFFREE(data->processing_list_1[i], processing_list_data);
            CC_IFFREE(data->processing_list_2[i], processing_list_data);
        }

        CC_IFFREE(data->g, GList**);
        CC_IFFREE(data->W, int*);
        CC_IFFREE(data->processing_list_1, processing_list_data*);
        CC_IFFREE(data->processing_list_2, processing_list_data*);
        CC_IFFREE(data, local_search_data);
    }

    return data;
}

void local_search_data_free(local_search_data** data) {
    if (*data != (local_search_data*)NULL) {
        printf("old iterations %d\n", (*data)->iterations);
        int nmachines = (*data)->nmachines;
        int njobs = (*data)->njobs;

        for (int i = 0; i < nmachines; ++i) {
            for (int j = 0; j < njobs; ++j) {
                g_list_free_full((*data)->g[i][j], destroy_slope_t);
            }

            CC_IFFREE((*data)->g[i], GList*);
            CC_IFFREE((*data)->W[i], int);
            CC_IFFREE((*data)->processing_list_1[i], processing_list_data);
            CC_IFFREE((*data)->processing_list_2[i], processing_list_data);
            CC_IFFREE((*data)->processing_list_1[i], processing_list_data);
            CC_IFFREE((*data)->processing_list_2[i], processing_list_data);
        }

        CC_IFFREE((*data)->g, GList**);
        CC_IFFREE((*data)->W, int*);
        CC_IFFREE((*data)->processing_list_1, processing_list_data*);
        CC_IFFREE((*data)->processing_list_2, processing_list_data*);
        CC_IFFREE((*data)->processing_list_1, processing_list_data*);
        CC_IFFREE((*data)->processing_list_2, processing_list_data*);
        CC_IFFREE(*data, local_search_data);
    }
}

int local_search_create_W(Solution* sol, local_search_data* data) {
    int  val = 0;
    int  nmachines;
    Job* tmp;

    if (sol == NULL || data == NULL || sol->nb_machines != data->nmachines) {
        val = 1;
        printf("Not compatible data structures\n");
        return val;
    }

    nmachines = data->nmachines;

    for (int i = 0; i < nmachines; ++i) {
        if (sol->part[i].used == 0) {
            continue;
        }

        tmp = (Job*)g_ptr_array_index(sol->part[i].machine, 0);
        data->W[i][0] = value_Fj(sol->c[tmp->job], tmp);

        for (unsigned j = 1; j < sol->part[i].machine->len; j++) {
            tmp = (Job*)g_ptr_array_index(sol->part[i].machine, j);
            data->W[i][j] = data->W[i][j - 1] + value_Fj(sol->c[tmp->job], tmp);
        }
    }

    return val;
}

static int local_search_create_processing_list(Solution*          sol,
                                               local_search_data* data,
                                               int                l) {
    int val = 0;

    for (int i = 0; i < data->nmachines; ++i) {
        int        njobs = sol->part[i].machine->len;
        int        C = 0;
        GPtrArray* machine = sol->part[i].machine;

        for (int j = 0; j < l; ++j) {
            Job* j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = 0; j < njobs - l; ++j) {
            Job* j1 = (Job*)g_ptr_array_index(machine, j);
            Job* j2 = (Job*)g_ptr_array_index(machine, j + l);
            data->processing_list_1[i][j].pos = j;
            data->processing_list_1[i][j].p = C;
            C = C - j1->processing_time + j2->processing_time;
        }

        qsort(data->processing_list_1[i], njobs - l,
              sizeof(processing_list_data), compare_process_list);
    }

    return val;
}

static int local_search_create_processing_list_2(Solution*          sol,
                                                 local_search_data* data,
                                                 int                l) {
    int val = 0;

    for (int i = 0; i < data->nmachines; ++i) {
        int        C = 0;
        int        njobs = sol->part[i].machine->len;
        GPtrArray* machine = sol->part[i].machine;

        for (int j = njobs - l; j < njobs; ++j) {
            Job* j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = njobs - l; j > 0; --j) {
            Job* j1 = (Job*)g_ptr_array_index(machine, j - 1);
            Job* j2 = (Job*)g_ptr_array_index(sol->part[i].machine, j + l - 1);
            data->processing_list_2[i][njobs - l - j].pos = j;
            data->processing_list_2[i][njobs - l - j].p = C;
            C = C + j1->processing_time - j2->processing_time;
        }

        qsort(data->processing_list_2[i], njobs - l,
              sizeof(processing_list_data), compare_process_list_b);
    }

    return val;
}

static int local_search_create_processing_list_swap(Solution*          sol,
                                                    local_search_data* data,
                                                    int                l1,
                                                    int                l2) {
    int val = 0;

    for (int i = 0; i < data->nmachines; ++i) {
        GPtrArray* machine = sol->part[i].machine;
        Job *      j1, *j2;
        int        njobs = sol->part[i].machine->len;
        int        C = 0;

        for (int j = l1; j < l1 + l2; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = l1; j < njobs - l2; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            j2 = (Job*)g_ptr_array_index(machine, j + l2);
            data->processing_list_2[i][j - l1].pos = j;
            data->processing_list_2[i][j - l1].p = C;
            C = C - j1->processing_time + j2->processing_time;
        }

        data->processing_list_2[i][njobs - l2 - l1].pos = njobs - l2;
        data->processing_list_2[i][njobs - l2 - l1].p = C;
        qsort(data->processing_list_2[i], njobs - l1 - l2 + 1,
              sizeof(processing_list_data), compare_process_list_b);
        C = 0;

        for (int j = 0; j < l1; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = 0; j < njobs - l1 - l2; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            j2 = (Job*)g_ptr_array_index(machine, j + l1);
            data->processing_list_1[i][j].pos = j;
            data->processing_list_1[i][j].p = C;
            C = C - j1->processing_time + j2->processing_time;
        }

        data->processing_list_1[i][njobs - l1 - l2].pos = njobs - l1 - l2;
        data->processing_list_1[i][njobs - l1 - l2].p = C;
        qsort(data->processing_list_1[i], njobs - l1 - l2 + 1,
              sizeof(processing_list_data), compare_process_list);
    }

    return val;
}

static int local_search_create_processing_list_insertion_inter(
    Solution*          sol,
    local_search_data* data,
    int                l) {
    int val = 0;

    for (int i = 0; i < data->nmachines; ++i) {
        int        njobs = sol->part[i].machine->len;
        int        C = 0;
        GPtrArray* machine = sol->part[i].machine;

        for (int j = 0; j < l; ++j) {
            Job* j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = 0; j < njobs - l; ++j) {
            Job* j1 = (Job*)g_ptr_array_index(machine, j);
            Job* j2 = (Job*)g_ptr_array_index(machine, j + l);
            data->processing_list_1[i][j].pos = j;
            data->processing_list_1[i][j].p = C;
            C = C - j1->processing_time + j2->processing_time;
        }

        data->processing_list_1[i][njobs - l].pos = njobs - l;
        data->processing_list_1[i][njobs - l].p = C;
        qsort(data->processing_list_1[i], njobs - l + 1,
              sizeof(processing_list_data), compare_process_list_b);
    }

    return val;
}

static int local_search_create_processing_list_swap_inter(
    Solution*          sol,
    local_search_data* data,
    int                l1,
    int                l2) {
    int val = 0;

    for (int i = 0; i < data->nmachines; ++i) {
        int        njobs = sol->part[i].machine->len;
        int        C = 0;
        Job *      j1, *j2;
        GPtrArray* machine = sol->part[i].machine;

        if (njobs - l1 < 0 || njobs - l2 < 0) {
            continue;
        }

        for (int j = 0; j < l1; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = 0; j < njobs - l1; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            j2 = (Job*)g_ptr_array_index(machine, j + l1);
            data->processing_list_1[i][j].pos = j;
            data->processing_list_1[i][j].p = C;
            C = C - j1->processing_time + j2->processing_time;
        }

        data->processing_list_1[i][njobs - l1].pos = njobs - l1;
        data->processing_list_1[i][njobs - l1].p = C;
        qsort(data->processing_list_1[i], njobs - l1 + 1,
              sizeof(processing_list_data), compare_process_list_b);
        C = 0;

        for (int j = 0; j < l2; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            C += j1->processing_time;
        }

        for (int j = 0; j < njobs - l2; ++j) {
            j1 = (Job*)g_ptr_array_index(machine, j);
            j2 = (Job*)g_ptr_array_index(machine, j + l2);
            data->processing_list_2[i][j].pos = j;
            data->processing_list_2[i][j].p = C;
            C = C - j1->processing_time + j2->processing_time;
        }

        data->processing_list_2[i][njobs - l2].pos = njobs - l2;
        data->processing_list_2[i][njobs - l2].p = C;
        qsort(data->processing_list_2[i], njobs - l2 + 1,
              sizeof(processing_list_data), compare_process_list_b);
    }

    return val;
}

int local_search_create_g(Solution* sol, local_search_data* data) {
    int  val = 0;
    int  nmachines = sol->nb_machines;
    int  njobs = sol->nb_jobs;
    Job* tmp;
    int  t1, t2;
    int  tw;
    int  w;

    for (int i = 0; i < nmachines; ++i) {
        if (sol->part[i].used == 0) {
            continue;
        } else {
            for (int j = 0; j < njobs; ++j) {
                g_list_free_full(data->g[i][j], destroy_slope_t);
                data->g[i][j] = (GList*)NULL;
            }

            sol->part[i].used = 0;
        }

        int n_k = sol->part[i].machine->len;
        int P = 0;

        for (int j = 0; j < n_k; ++j) {
            GPtrArray* lateness_sort = g_ptr_array_new();

            for (int k = j; k < n_k; ++k) {
                g_ptr_array_add(lateness_sort,
                                g_ptr_array_index(sol->part[i].machine, k));
            }

            g_ptr_array_sort_with_data(lateness_sort,
                                       local_search_compare_lateness, sol->c);
            tw = 0;
            w = 0;
            t1 = 0;
            unsigned k = 0;
            int      move;
            move = 1;
            tmp = (Job*)g_ptr_array_index(lateness_sort, 0);

            for (; k < lateness_sort->len && move;) {
                move = tmp->weight * (sol->c[tmp->job] - P - tmp->due_time) > 0;

                if (move) {
                    tw += tmp->weight * (sol->c[tmp->job] - P - tmp->due_time);
                    w += tmp->weight;
                    k++;
                    tmp = (Job*)g_ptr_array_index(lateness_sort, k);
                }
            }

            t2 = tmp->due_time - sol->c[tmp->job] + P;
            local_search_add_slope_t(data, t1, t2, tw, w, i, j);

            for (unsigned l = k; l < lateness_sort->len;) {
                tw = tw + w * (t2 - t1);
                t1 = t2;
                move = 1;
                tmp = (Job*)g_ptr_array_index(lateness_sort, l);

                while (move) {
                    w += tmp->weight;
                    l++;

                    if (l == lateness_sort->len) {
                        move = 0;
                        t2 = INT_MAX;
                    } else {
                        tmp = (Job*)g_ptr_array_index(lateness_sort, l);
                        t2 = tmp->due_time - sol->c[tmp->job] + P;
                        move = (t1 == t2);
                    }
                }

                local_search_add_slope_t(data, t1, t2, tw, w, i, j);
            }

            tmp = (Job*)g_ptr_array_index(sol->part[i].machine, j);
            P += tmp->processing_time;
            g_ptr_array_free(lateness_sort, TRUE);
        }
    }

    return val;
}

static void local_search_update_insertion(Solution*        sol,
                                          int              i_best,
                                          int              j_best,
                                          int              k_best,
                                          int              l,
                                          MAYBE_UNUSED int improvement) {
    Job* tmp;
#ifndef NDEBUG
    int old = sol->tw;
#endif
    sol->tw -= sol->part[k_best].tw;
    sol->part[k_best].c = 0;
    sol->part[k_best].tw = 0;

    for (int i = 0; i < l; ++i) {
        tmp = (Job*)g_ptr_array_index(sol->part[k_best].machine, i_best);
        g_ptr_array_remove_index(sol->part[k_best].machine, i_best);
        g_ptr_array_insert(sol->part[k_best].machine, j_best, tmp);
    }

    for (unsigned i = 0; i < sol->part[k_best].machine->len; ++i) {
        tmp = (Job*)g_ptr_array_index(sol->part[k_best].machine, i);
        sol->part[k_best].c += tmp->processing_time;
        sol->c[tmp->job] = sol->part[k_best].c;
        sol->part[k_best].tw += value_Fj(sol->c[tmp->job], tmp);
    }

    sol->tw += sol->part[k_best].tw;
    sol->part[k_best].used = 1;
    assert(old - sol->tw == improvement);
}

static void local_search_update_insertion_inter(Solution*        sol,
                                                int              i_best,
                                                int              j_best,
                                                int              k_best,
                                                int              kk_best,
                                                int              l,
                                                MAYBE_UNUSED int improvement) {
    Job* tmp;
#ifndef NDEBUG
    int old = sol->tw;
#endif
    sol->tw -= sol->part[k_best].tw;
    sol->tw -= sol->part[kk_best].tw;
    sol->part[k_best].c = 0;
    sol->part[k_best].tw = 0;
    sol->part[kk_best].c = 0;
    sol->part[kk_best].tw = 0;

    for (int i = 0; i < l; ++i) {
        tmp = (Job*)g_ptr_array_index(sol->part[k_best].machine, i_best);
        g_ptr_array_remove_index(sol->part[k_best].machine, i_best);
        g_ptr_array_insert(sol->part[kk_best].machine, j_best + i, tmp);
    }

    for (unsigned i = 0; i < sol->part[k_best].machine->len; ++i) {
        tmp = (Job*)g_ptr_array_index(sol->part[k_best].machine, i);
        sol->part[k_best].c += tmp->processing_time;
        sol->c[tmp->job] = sol->part[k_best].c;
        sol->part[k_best].tw +=
            tmp->weight * CC_MAX(0, sol->c[tmp->job] - tmp->due_time);
    }

    for (unsigned i = 0; i < sol->part[kk_best].machine->len; ++i) {
        tmp = (Job*)g_ptr_array_index(sol->part[kk_best].machine, i);
        sol->part[kk_best].c += tmp->processing_time;
        sol->c[tmp->job] = sol->part[kk_best].c;
        sol->part[kk_best].tw += value_Fj(sol->c[tmp->job], tmp);
    }

    sol->tw += sol->part[k_best].tw + sol->part[kk_best].tw;
    assert(old - sol->tw == improvement);
    sol->part[k_best].used = 1;
    sol->part[kk_best].used = 1;
}

static void local_search_update_swap(Solution*        sol,
                                     int              i_best,
                                     int              j_best,
                                     int              k_best,
                                     int              l1,
                                     int              l2,
                                     MAYBE_UNUSED int improvement) {
    Job*      tmp;
    gpointer  swap;
    PartList* part = sol->part + k_best;
#ifndef NDEBUG
    int old = sol->tw;
#endif
    sol->tw -= part->tw;
    part->c = 0;
    part->tw = 0;

    if (l1 == l2) {
        for (int i = 0; i < l1; ++i) {
            CC_SWAP(g_ptr_array_index(part->machine, i_best + i),
                    g_ptr_array_index(part->machine, j_best + i), swap);
        }
    } else if (l1 < l2) {
        for (int i = 0; i < l1; ++i) {
            CC_SWAP(g_ptr_array_index(part->machine, i_best + i),
                    g_ptr_array_index(part->machine, j_best + i), swap);
        }

        for (int i = l1; i < l2; ++i) {
            tmp = (Job*)g_ptr_array_index(part->machine, j_best + i);
            g_ptr_array_remove_index(part->machine, j_best + i);
            g_ptr_array_insert(part->machine, i_best + i, tmp);
        }
    } else {
        for (int i = 0; i < l2; ++i) {
            CC_SWAP(g_ptr_array_index(part->machine, i_best + i),
                    g_ptr_array_index(part->machine, j_best + i), swap);
        }

        for (int i = l2; i < l1; ++i) {
            tmp = (Job*)g_ptr_array_index(part->machine, i_best + i);
            g_ptr_array_remove_index(part->machine, i_best + i);
            g_ptr_array_insert(part->machine, j_best + i - 1, tmp);
        }
    }

    for (unsigned i = 0; i < part->machine->len; ++i) {
        tmp = (Job*)g_ptr_array_index(sol->part[k_best].machine, i);
        part->c += tmp->processing_time;
        sol->c[tmp->job] = part->c;
        part->tw += value_Fj(sol->c[tmp->job], tmp);
    }

    sol->tw += part->tw;
    assert(old - sol->tw == improvement);
    part->used = 1;
}

static void local_search_update_inter_swap(Solution*        sol,
                                           int              i_best,
                                           int              j_best,
                                           int              k_best,
                                           int              kk_best,
                                           int              l1,
                                           int              l2,
                                           MAYBE_UNUSED int improvement) {
    Job*      tmp;
    gpointer  swap;
    PartList* part1 = sol->part + k_best;
    PartList* part2 = sol->part + kk_best;
#ifndef NDEBUG
    int old = sol->tw;
#endif
    sol->tw -= part1->tw + part2->tw;
    part1->c = 0;
    part1->tw = 0;
    part2->c = 0;
    part2->tw = 0;

    if (l1 == l2) {
        for (int i = 0; i < l1; ++i) {
            CC_SWAP(g_ptr_array_index(part1->machine, i_best + i),
                    g_ptr_array_index(part2->machine, j_best + i), swap);
        }
    } else if (l1 < l2) {
        for (int i = 0; i < l1; ++i) {
            CC_SWAP(g_ptr_array_index(part1->machine, i_best + i),
                    g_ptr_array_index(part2->machine, j_best + i), swap);
        }

        for (int i = l1; i < l2; ++i) {
            tmp = (Job*)g_ptr_array_index(part2->machine, j_best + l1);
            g_ptr_array_remove_index(part2->machine, j_best + l1);
            g_ptr_array_insert(part1->machine, i_best + i, tmp);
        }
    } else {
        for (int i = 0; i < l2; ++i) {
            CC_SWAP(g_ptr_array_index(part1->machine, i_best + i),
                    g_ptr_array_index(part1->machine, j_best + i), swap);
        }

        for (int i = l2; i < l1; ++i) {
            tmp = (Job*)g_ptr_array_index(part1->machine, i_best + l2);
            g_ptr_array_remove_index(part1->machine, i_best + l2);
            g_ptr_array_insert(part2->machine, j_best + i, tmp);
        }
    }

    for (unsigned i = 0; i < part1->machine->len; ++i) {
        tmp = (Job*)g_ptr_array_index(part1->machine, i);
        part1->c += tmp->processing_time;
        sol->c[tmp->job] = part1->c;
        part1->tw += value_Fj(sol->c[tmp->job], tmp);
    }

    for (unsigned i = 0; i < part2->machine->len; ++i) {
        tmp = (Job*)g_ptr_array_index(part2->machine, i);
        part2->c += tmp->processing_time;
        sol->c[tmp->job] = part2->c;
        part2->tw += value_Fj(sol->c[tmp->job], tmp);
    }

    sol->tw += part1->tw + part2->tw;
    assert(old - sol->tw == improvement);
    part1->used = 1;
    part2->used = 1;
}

void local_search_forward_insertion(Solution*          sol,
                                    local_search_data* data,
                                    int                l) {
    int    pos, p, c, tmp;
    int    update;
    Job*   tmp_j;
    int    max;
    int    i_best = -1, j_best = -1, k_best = -1;
    double runningtime = CCutil_zeit();

    update = 0;
    max = 0;
    data->iterations++;
    local_search_create_processing_list(sol, data, l);

    for (int k = 0; k < sol->nb_machines; ++k) {
        /** compute g */
        int njobs = sol->part[k].machine->len;

        for (int i = 0; i < njobs - l; ++i) {
            pos = data->processing_list_1[k][i].pos;
            p = data->processing_list_1[k][i].p;
            GList* it = data->g[k][pos];

            for (int j = pos + l; j < njobs; ++j) {
                tmp_j = (Job*)g_ptr_array_index(sol->part[k].machine, j);
                c = sol->c[tmp_j->job] - p;
                compute_it(&it, c);
                g[pos][j] = compute_g(&it, c);
            }
        }

        /** compute h */
        for (int i = 0; i < njobs - l; ++i) {
            GList* it = data->g[k][i + l];

            for (int j = i + l; j < njobs; ++j) {
                tmp_j = (Job*)g_ptr_array_index(sol->part[k].machine, j);
                c = sol->c[tmp_j->job];
                compute_it(&it, c);
                h[i][j] = compute_g(&it, sol->c[tmp_j->job]);
            }
        }

        /** compute gg */
        for (int i = 0; i < njobs - l; ++i) {
            GList* it = data->g[k][i + l];

            if (i != 0) {
                tmp_j = (Job*)g_ptr_array_index(sol->part[k].machine, i - 1);
                c = sol->c[tmp_j->job];
            } else {
                c = 0;
            }

            compute_it(&it, c);
            gg[i] = compute_g(&it, c);
        }

        /** compute hh */
        for (int j = l; j < njobs - 1; ++j) {
            iterators[j] = data->g[k][j + 1];
        }

        iterators[njobs - 1] = (GList*)NULL;

        for (int i = 0; i < njobs - l; ++i) {
            pos = data->processing_list_1[k][i].pos;
            p = data->processing_list_1[k][i].p;

            for (int j = pos + l; j < njobs; ++j) {
                if (iterators[j] != (GList*)NULL) {
                    tmp_j = (Job*)g_ptr_array_index(sol->part[k].machine, j);
                    c = sol->c[tmp_j->job] - p;
                    compute_it(&iterators[j], c);
                    hh[pos][j] = compute_g(&iterators[j], c);
                } else {
                    hh[pos][j] = 0;
                }
            }
        }

        for (int i = 0; i < njobs - l; ++i) {
            for (int j = i + l; j < njobs; ++j) {
                tmp = 0;

                if (i != 0) {
                    tmp += data->W[k][i - 1];
                }

                tmp += g[i][j] - h[i][j];
                tmp += gg[i] - hh[i][j];
                tmp += data->W[k][njobs - 1] - data->W[k][j];

                if (sol->part[k].tw - tmp > max) {
                    max = sol->part[k].tw - tmp;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    update = 1;
                }
            }
        }
    }

    /** update to best improvement */
    if (update) {
        if (dbg_lvl()) {
            solution_print(sol);
        }

        local_search_update_insertion(sol, i_best, j_best, k_best, l, max);
        local_search_create_W(sol, data);
        local_search_create_g(sol, data);

        if (dbg_lvl()) {
            solution_print(sol);
        }

        data->updated = 1;
    } else {
        data->updated = 0;
    }

    if (dbg_lvl() > 0) {
        printf(
            "forward insertion with l = %d, running time = %f and improvement "
            "%d\n",
            l, CCutil_zeit() - runningtime, max);
        print_line();
    }
}

void local_search_backward_insertion(Solution*          sol,
                                     local_search_data* data,
                                     int                l) {
    int    c;
    int    pos;
    int    t;
    int    update;
    GList* it;
    Job*   tmp_j;
    int    max;
    int    i_best = -1, j_best = -1, k_best = -1;
    double runningtime = CCutil_zeit();

    update = 0;
    max = 0;
    local_search_create_processing_list_2(sol, data, l);

    data->iterations++;
    for (int k = 0; k < sol->nb_machines; ++k) {
        int        njobs = sol->part[k].machine->len;
        int        p;
        GPtrArray* machine = sol->part[k].machine;

        /** compute g */
        for (int i = njobs - l; i > 0; --i) {
            it = data->g[k][i];

            c = 0;
            for (int j = 1; j < i; ++j) {
                tmp_j = (Job*)g_ptr_array_index(machine, j - 1);
                c = sol->c[tmp_j->job];

                compute_it(&it, c);
                g[i][j] = compute_g(&it, c);
            }
        }

        /** compute h */
        for (int j = 0; j < njobs - l; ++j) {
            h[njobs - l][j] = 0;
        }

        p = 0;

        for (int i = njobs - l - 1; i < njobs - 1; ++i) {
            tmp_j = (Job*)g_ptr_array_index(machine, i);
            p += tmp_j->processing_time;
        }

        for (int i = njobs - l - 1; i > 0; --i) {
            it = data->g[k][i + l];

            c = p;
            for (int j = 1; j < i; ++j) {
                tmp_j = (Job*)g_ptr_array_index(machine, j - 1);
                c = p + sol->c[tmp_j->job];

                compute_it(&it, c);
                h[i][j] = compute_g(&it, c);
            }

            tmp_j = (Job*)g_ptr_array_index(machine, i + l - 1);
            p -= tmp_j->processing_time;
            tmp_j = (Job*)g_ptr_array_index(machine, i - 1);
            p += tmp_j->processing_time;
        }

        /** compute gg */
        for (int i = njobs - l; i > 0; --i) {
            it = data->g[k][i];
            tmp_j = (Job*)g_ptr_array_index(machine, i + l - 1);
            c = sol->c[tmp_j->job];
            compute_it(&it, c);
            gg[i] = compute_g(&it, c);
        }

        /** compute hh */
        for (int j = 0; j < njobs - l; ++j) {
            iterators[j] = data->g[k][j];
        }

        for (int i = njobs - l - 1; i > 0; --i) {
            p = data->processing_list_2[k][njobs - l - i].p;
            pos = data->processing_list_2[k][njobs - l - i].pos;

            c = p;
            for (int j = 1; j < pos; ++j) {
                tmp_j = (Job*)g_ptr_array_index(machine, j - 1);
                c = p + sol->c[tmp_j->job];

                compute_it(&iterators[j], c);
                hh[pos][j] = compute_g(&iterators[j], c);
            }
        }

        for (int i = njobs - l; i > 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                t = 0;

                // if (j != 0) {
                t += j != 0 ? data->W[k][j - 1] : 0;
                // }

                t += g[i][j] - h[i][j];
                t += hh[i][j] - gg[i];
                t += data->W[k][njobs - 1] - data->W[k][i + l - 1];

                if (sol->part[k].tw - t > max) {
                    max = sol->part[k].tw - t;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    update = 1;
                }
            }
        }
    }

    /** update to best improvement */
    if (update) {
        if (dbg_lvl() > 0) {
            solution_print(sol);
        }

        local_search_update_insertion(sol, i_best, j_best, k_best, l, max);
        local_search_create_W(sol, data);
        local_search_create_g(sol, data);

        if (dbg_lvl() > 0) {
            solution_print(sol);
        }

        data->updated = 1;
    } else {
        data->updated = 0;
    }

    if (dbg_lvl() > 0) {
        printf(
            "backward insertion with l = %d, running time = %f and improvement "
            "%d\n",
            l, CCutil_zeit() - runningtime, max);
        print_line();
    }
}

void local_search_swap_intra(Solution*          sol,
                             local_search_data* data,
                             int                l1,
                             int                l2) {
    int    pos, p, c, t;
    int    update;
    GList* it;
    Job*   tmp_j;
    int    max;
    int    i_best = -1, j_best = -1, k_best = -1;
    double runningtime = CCutil_zeit();

    update = 0;
    max = 0;
    local_search_create_processing_list_swap(sol, data, l1, l2);
    data->iterations++;
    for (int k = 0; k < sol->nb_machines; ++k) {
        int        njobs = sol->part[k].machine->len;
        GPtrArray* machine = sol->part[k].machine;

        /** compute g */
        for (int i = 0; i < njobs - l1 - l2 + 1; ++i) {
            pos = data->processing_list_1[k][i].pos;
            p = data->processing_list_1[k][i].p;
            it = data->g[k][pos];

            for (int j = pos + l1; j < njobs - l2 + 1; ++j) {
                tmp_j = (Job*)g_ptr_array_index(machine, j + l2 - 1);
                c = sol->c[tmp_j->job] - p;
                compute_it(&it, c);
                B2_1[pos][j] = compute_g(&it, c);
            }
        }

        /** compute h */
        for (int i = 0; i < njobs - l1 - l2 + 1; ++i) {
            it = (i + l1 >= njobs) ? NULL : data->g[k][i + l1];

            for (int j = i + l1; j < njobs - l2 + 1; ++j) {
                if (it == NULL) {
                    B2_2[i][j] = 0;
                } else {
                    tmp_j = (Job*)g_ptr_array_index(machine, j + l2 - 1);
                    c = sol->c[tmp_j->job];
                    compute_it(&it, c);
                    B2_2[i][j] = compute_g(&it, c);
                }
            }
        }

        /** compute B3_1 */
        for (int i = 0; i < njobs - l1 - l2 + 1; ++i) {
            iterators[i] = data->g[k][i + l1];
        }

        for (int j = l1; j < njobs - l2 + 1; ++j) {
            pos = data->processing_list_2[k][j - l1].pos;
            p = data->processing_list_2[k][j - l1].p;

            for (int i = 0; i < pos - l1 + 1; ++i) {
                if (i + l1 >= njobs) {
                    B3_1[i][pos] = 0;
                } else {
                    c = p;

                    if (i != 0) {
                        tmp_j = (Job*)g_ptr_array_index(machine, i - 1);
                        c += sol->c[tmp_j->job];
                    }

                    compute_it(&iterators[i], c);
                    B3_1[i][pos] = compute_g(&iterators[i], c);
                }
            }
        }

        /** compute B3_2 */
        for (int j = l1; j < njobs - l2 + 1; ++j) {
            iterators[j] = data->g[k][j];
        }

        for (int i = 0; i < njobs - l1 - l2 + 1; ++i) {
            pos = data->processing_list_1[k][i].pos;
            p = data->processing_list_1[k][i].p;

            for (int j = pos + l1; j < njobs - l2 + 1; ++j) {
                if (i + l1 >= njobs) {
                    B3_2[pos][j] = 0;
                } else {
                    tmp_j = (Job*)g_ptr_array_index(machine, j + l2 - 1);
                    c = sol->c[tmp_j->job] - p;
                    compute_it(&iterators[j], c);
                    B3_2[pos][j] = compute_g(&iterators[j], c);
                }
            }
        }

        /** compute B4_1 */
        for (int j = l1; j < njobs - l2 + 1; ++j) {
            it = data->g[k][j];

            for (int i = 0; i < j - l1 + 1; ++i) {
                c = 0;

                if (i != 0) {
                    c = sol->c[((Job*)g_ptr_array_index(machine, i - 1))->job];
                }

                compute_it(&it, c);
                B4_1[i][j] = compute_g(&it, c);
            }
        }

        /** compute B4_2 */
        for (int j = l1; j < njobs - l2 + 1; ++j) {
            pos = data->processing_list_2[k][j - l1].pos;
            p = data->processing_list_2[k][j - l1].p;

            if (pos + l2 == njobs) {
                it = NULL;
            } else {
                it = data->g[k][pos + l2];
            }

            for (int i = 0; i < pos - l1 + 1; ++i) {
                c = p;

                if (i != 0) {
                    c += sol->c[((Job*)g_ptr_array_index(machine, i - 1))->job];
                }

                compute_it(&it, c);

                if (it != NULL) {
                    B4_2[i][pos] = compute_g(&it, c);
                } else {
                    B4_2[i][pos] = 0;
                }
            }
        }

        for (int i = 0; i < njobs - l1 - l2 + 1; ++i) {
            for (int j = i + l1; j < njobs - l2 + 1; ++j) {
                t = 0;

                if (i != 0) {
                    t += data->W[k][i - 1];
                }

                t += B2_1[i][j] - B2_2[i][j];
                t += B3_1[i][j] - B3_2[i][j];
                t += B4_1[i][j] - B4_2[i][j];
                t += data->W[k][njobs - 1] - data->W[k][j + l2 - 1];

                if (sol->part[k].tw - t > max) {
                    max = sol->part[k].tw - t;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    update = 1;
                }
            }
        }
    }

    /** update to best improvement */
    if (update) {
        if (dbg_lvl()) {
            solution_print(sol);
        }

        local_search_update_swap(sol, i_best, j_best, k_best, l1, l2, max);
        local_search_create_W(sol, data);
        local_search_create_g(sol, data);

        if (dbg_lvl()) {
            solution_print(sol);
        }

        data->updated = 1;
    } else {
        data->updated = 0;
    }

    if (dbg_lvl() > 0) {
        printf(
            "intra swap with l1 = %d and l2 = %d, running time = %f and "
            "improvement %d on machine %d on places %d %d\n",
            l1, l2, CCutil_zeit() - runningtime, max, k_best, i_best, j_best);
        print_line();
    }
}

void local_search_insertion_inter(Solution*          sol,
                                  local_search_data* data,
                                  int                l) {
    int    pos, p, c, t;
    int    update;
    GList* it;
    Job*   tmp_j;
    int    max;
    int    i_best = -1, j_best = -1, k_best = -1, kk_best = -1;
    double runningtime = CCutil_zeit();

    update = 0;
    max = 0;
    local_search_create_processing_list_insertion_inter(sol, data, l);
    data->iterations++;
    for (int k1 = 0; k1 < sol->nb_machines; ++k1) {
        int        njobs1 = sol->part[k1].machine->len;
        GPtrArray* machine1 = sol->part[k1].machine;

        for (int k2 = 0; k2 < sol->nb_machines; ++k2) {
            if (k1 == k2) {
                continue;
            }

            int        njobs2 = sol->part[k2].machine->len;
            GPtrArray* machine2 = sol->part[k2].machine;

            /** compute B2_1 */
            for (int i = 0; i < njobs1 - l + 1; ++i) {
                it = data->g[k1][i];

                for (int j = 0; j < njobs2; ++j) {
                    c = 0;

                    if (j != 0) {
                        tmp_j = (Job*)g_ptr_array_index(machine2, j - 1);
                        c = sol->c[tmp_j->job];
                    }

                    compute_it(&it, c);
                    B2_1[i][j] = compute_g(&it, c);
                }
            }

            for (int i = 0; i < njobs1 - l + 1; ++i) {
                pos = data->processing_list_1[k1][i].pos;
                p = data->processing_list_1[k1][i].p;
                it = (pos + l >= njobs1) ? NULL : data->g[k1][pos + l];

                for (int j = 0; j < njobs2; ++j) {
                    if (it == NULL) {
                        B2_2[pos][j] = 0;
                    } else {
                        c = p;

                        if (j != 0) {
                            tmp_j = (Job*)g_ptr_array_index(machine2, j - 1);
                            c += sol->c[tmp_j->job];
                        }

                        compute_it(&it, c);
                        B2_2[pos][j] = compute_g(&it, c);
                    }
                }
            }

            for (int i = 0; i < njobs1 - l + 1; ++i) {
                if (i + l >= njobs1) {
                    B3_1_[i] = 0;
                } else {
                    it = data->g[k1][i];
                    c = 0;

                    if (i != 0) {
                        tmp_j = (Job*)g_ptr_array_index(machine1, i - 1);
                        c = sol->c[tmp_j->job];
                    }

                    compute_it(&it, c);
                    B3_1_[i] = compute_g(&it, c);
                }
            }

            for (int j = 0; j < njobs2 - 1; ++j) {
                iterators[j] = data->g[k2][j];
            }

            for (int j = 0; j < njobs2 - 1; ++j) {
                it = data->g[k2][j];

                for (int i = 0; i < njobs1 - l + 1; ++i) {
                    pos = data->processing_list_1[k1][i].pos;
                    p = data->processing_list_1[k1][i].p;
                    c = p;

                    if (j != 0) {
                        tmp_j = (Job*)g_ptr_array_index(machine2, j - 1);
                        c += sol->c[tmp_j->job];
                    }

                    compute_it(&it, c);
                    B5_1[pos][j] = compute_g(&it, c);
                }
            }

            for (int i = 0; i < njobs1 - l + 1; ++i) {
                for (int j = 0; j < njobs2 - 1; ++j) {
                    t = 0;

                    if (i != 0) {
                        t += data->W[k1][i - 1];
                    }

                    t += B2_1[i][j] - B2_2[i][j];
                    t += B3_1_[i];
                    t += B5_1[i][j];

                    if (j != 0) {
                        t += data->W[k2][j - 1];
                    }

                    if (sol->part[k1].tw + sol->part[k2].tw - t > max) {
                        max = sol->part[k1].tw + sol->part[k2].tw - t;
                        i_best = i;
                        j_best = j;
                        k_best = k1;
                        kk_best = k2;
                        update = 1;
                    }
                }
            }
        }
    }

    /** update to best improvement */
    if (update) {
        if (dbg_lvl()) {
            solution_print(sol);
        }

        local_search_update_insertion_inter(sol, i_best, j_best, k_best,
                                            kk_best, l, max);
        local_search_create_W(sol, data);
        local_search_create_g(sol, data);
        data->updated = 1;

        if (dbg_lvl()) {
            solution_print(sol);
        }
    } else {
        data->updated = 0;
    }

    if (dbg_lvl() > 0) {
        printf(
            "inter insertion with l = %d, running time = %f and improvement %d "
            "on machines %d and %d on places %d %d\n",
            l, CCutil_zeit() - runningtime, max, k_best, kk_best, i_best,
            j_best);
        print_line();
    }
}

void local_search_swap_inter(Solution*          sol,
                             local_search_data* data,
                             int                l1,
                             int                l2) {
    int    pos, p, c, t;
    int    update;
    GList* it;
    Job*   tmp_j;
    int    max;
    int    i_best = -1, j_best = -1, k_best = -1, kk_best = -1;
    double runningtime = CCutil_zeit();

    update = 0;
    max = 0;
    local_search_create_processing_list_swap_inter(sol, data, l1, l2);
    data->iterations++;
    for (int k1 = 0; k1 < sol->nb_machines; ++k1) {
        int        njobs1 = sol->part[k1].machine->len;
        GPtrArray* machine1 = sol->part[k1].machine;

        for (int k2 = 0; k2 < sol->nb_machines; ++k2) {
            int        njobs2 = sol->part[k2].machine->len;
            GPtrArray* machine2 = sol->part[k2].machine;

            if (k1 == k2 || njobs1 - l1 < 0 || njobs2 - l2 < 0) {
                continue;
            }

            /** compute B2_1 */
            for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                it = data->g[k1][i];

                for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                    c = 0;

                    if (j != 0) {
                        tmp_j = (Job*)g_ptr_array_index(machine2, j - 1);
                        c = sol->c[tmp_j->job];
                    }

                    compute_it(&it, c);
                    B2_1[i][j] = compute_g(&it, c);
                }
            }

            for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                pos = data->processing_list_1[k1][i].pos;
                p = data->processing_list_1[k1][i].p;
                it = (pos + l1 >= njobs1) ? NULL : data->g[k1][pos + l1];

                for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                    if (it == NULL) {
                        B2_2[pos][j] = 0;
                    } else {
                        c = p;

                        if (j != 0) {
                            tmp_j = (Job*)g_ptr_array_index(machine2, j - 1);
                            c += sol->c[tmp_j->job];
                        }

                        compute_it(&it, c);
                        B2_2[pos][j] = compute_g(&it, c);
                    }
                }
            }

            for (int i = 0; i < njobs1 - l1; ++i) {
                iterators[i] = data->g[k1][i + l1];
            }

            for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                pos = data->processing_list_2[k2][j].pos;
                p = data->processing_list_2[k2][j].p;

                for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                    if (i + l1 >= njobs1) {
                        B3_1[i][pos] = 0;
                    } else {
                        c = p;

                        if (i != 0) {
                            tmp_j = (Job*)g_ptr_array_index(machine1, i - 1);
                            c += sol->c[tmp_j->job];
                        }

                        compute_it(&iterators[i], c);
                        B3_1[i][pos] = compute_g(&iterators[i], c);
                    }
                }
            }

            for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                it = data->g[k2][j];

                for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                    c = 0;

                    if (i != 0) {
                        tmp_j = (Job*)g_ptr_array_index(machine1, i - 1);
                        c += sol->c[tmp_j->job];
                    }

                    compute_it(&it, c);
                    B5_1[i][j] = compute_g(&it, c);
                }
            }

            for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                pos = data->processing_list_2[k2][j].pos;
                p = data->processing_list_2[k2][j].p;
                it = (pos + l2 >= njobs2) ? NULL : data->g[k2][pos + l2];

                for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                    if (it == NULL) {
                        B5_2[i][pos] = 0;
                    } else {
                        c = p;

                        if (i != 0) {
                            tmp_j = (Job*)g_ptr_array_index(machine1, i - 1);
                            c += sol->c[tmp_j->job];
                        }

                        compute_it(&it, c);
                        B5_2[i][pos] = compute_g(&it, c);
                    }
                }
            }

            for (int j = 0; j < njobs2 - l2; ++j) {
                iterators[j] = data->g[k2][j + l2];
            }

            for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                pos = data->processing_list_1[k1][i].pos;
                p = data->processing_list_1[k1][i].p;

                for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                    if (j + l2 >= njobs2) {
                        B6_1[pos][j] = 0;
                    } else {
                        c = p;

                        if (j != 0) {
                            tmp_j = (Job*)g_ptr_array_index(machine2, j - 1);
                            c += sol->c[tmp_j->job];
                        }

                        compute_it(&iterators[j], c);
                        B6_1[pos][j] = compute_g(&iterators[j], c);
                    }
                }
            }

            for (int i = 0; i < njobs1 - l1 + 1; ++i) {
                for (int j = 0; j < njobs2 - l2 + 1; ++j) {
                    t = 0;

                    if (i != 0) {
                        t += data->W[k1][i - 1];
                    }

                    t += B2_1[i][j] - B2_2[i][j];
                    t += B3_1[i][j];
                    t += B5_1[i][j] - B5_2[i][j];
                    t += B6_1[i][j];

                    if (j != 0) {
                        t += data->W[k2][j - 1];
                    }

                    if (sol->part[k1].tw + sol->part[k2].tw - t > max) {
                        max = sol->part[k1].tw + sol->part[k2].tw - t;
                        i_best = i;
                        j_best = j;
                        k_best = k1;
                        kk_best = k2;
                        update = 1;
                    }
                }
            }
        }
    }

    /** update to best improvement */
    if (update) {
        if (dbg_lvl()) {
            solution_print(sol);
        }

        local_search_update_inter_swap(sol, i_best, j_best, k_best, kk_best, l1,
                                       l2, max);
        local_search_create_W(sol, data);
        local_search_create_g(sol, data);
        data->updated = 1;

        if (dbg_lvl()) {
            solution_print(sol);
        }
    } else {
        data->updated = 0;
    }

    if (dbg_lvl() > 0) {
        printf(
            "inter insertion with l1 = %d and l2 = %d, running time = %f and "
            "improvement %d on machines %d and %d on places %d %d\n",
            l1, l2, CCutil_zeit() - runningtime, max, k_best, kk_best, i_best,
            j_best);
        print_line();
    }
}
