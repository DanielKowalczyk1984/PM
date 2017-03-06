#include <assert.h>
#include "wct.h"

void destroy_slope_t(gpointer data){
    slope_t *tmp = (slope_t *) data;
    CC_IFFREE(tmp, slope_t);
}

int compute_g(GList **it, int t) {
    slope_t *x = (slope_t *) (*it)->data;
    return x->c + x->alpha * (t - x->b1);
}

void compute_it(GList **it, int c){
    if((*it) != NULL) {
        slope_t * tmp = (slope_t *) (*it)->data;
        int move = !(tmp->b1 <= c &&
                 tmp->b2 >= c);

        while (move) {
            *it = (*it)->next;

            if (*it != (GList *) NULL) {
                tmp = (slope_t *) (*it)->data;
                move = !(tmp->b1 <= c &&
                         tmp->b2 > c);
            } else {
                move = 0;
            }
        }
    }
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

int compare_process_list_b(gconstpointer a, gconstpointer b) {
    processing_list_data *x = (processing_list_data *) a;
    processing_list_data *y = (processing_list_data *) b;

    if (x->p > y->p) {
        return 1;
    } else if (x->p < y->p) {
        return -1;
    }

    return 0;

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

static void local_search_add_slope_t(local_search_data *data, int b1, int b2,
                                     int c, int alpha, int i, int j) {
    slope_t *tmp = CC_SAFE_MALLOC(1, slope_t);
    tmp->alpha = alpha;
    tmp->c = c;
    tmp->b1 = b1;
    tmp->b2 = b2;
    data->g[i][j] = g_list_append(data->g[i][j],tmp);
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
    data->processing_list_f = CC_SAFE_MALLOC(nmachines, processing_list_data *);
    data->processing_list_b = CC_SAFE_MALLOC(nmachines, processing_list_data *);
    data->njobs = CC_SAFE_MALLOC(nmachines, int);

    for (i = 0; i < nmachines; ++i) {
        data->njobs[i] = sol->part[i].machine->len;
        data->W[i] = CC_SAFE_MALLOC(data->njobs[i], int);
        data->g[i] = CC_SAFE_MALLOC(data->njobs[i], GList *);
        data->processing_list_f[i] = (processing_list_data *) NULL;
        data->processing_list_b[i] = (processing_list_data *) NULL;

        for (j = 0; j < sol->part[i].machine->len; ++j) {
            data->g[i][j] = (GList *) NULL;
        }
    }


CLEAN:

    if (val) {
        for (i = 0; i < nmachines; ++i) {
            for (j = 0; j < sol->part[i].machine->len; ++j) {
                g_list_free_full(data->g[i][j],destroy_slope_t);
            }

            CC_IFFREE(data->g[i], GList *);
            CC_IFFREE(data->W[i], int);

            if (data->processing_list_f[i] != NULL) {
                CC_IFFREE(data->processing_list_f[i], processing_list_data);
                CC_IFFREE(data->processing_list_b[i], processing_list_data);
            }
        }

        CC_IFFREE(data->g, GList **);
        CC_IFFREE(data->W, int *);
        CC_IFFREE(data->processing_list_f, processing_list_data *);
        CC_IFFREE(data->processing_list_b, processing_list_data *);
        CC_IFFREE(data->njobs, int);
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
                g_list_free_full(data->g[i][j], destroy_slope_t);
            }

            CC_IFFREE(data->g[i], GList *);
            CC_IFFREE(data->W[i], int);

            if (data->processing_list_f[i] != NULL) {
                CC_IFFREE(data->processing_list_f[i], processing_list_data);
                CC_IFFREE(data->processing_list_b[i], processing_list_data);
            }
        }

        CC_IFFREE(data->g, GList **);
        CC_IFFREE(data->W, int *);
        CC_IFFREE(data->processing_list_f, processing_list_data *);
        CC_IFFREE(data->processing_list_b, processing_list_data *);
        CC_IFFREE(data->njobs, int);
        CC_IFFREE(data, local_search_data);
    }
}



int local_search_create_W(solution *sol, local_search_data *data) {
    int val = 0;
    int nmachines;
    Job *tmp;

    if (sol == NULL || data == NULL || sol->nmachines != data->nmachines) {
        val = 1;
        printf("Not compatible data structures\n");
        return val;
    }

    nmachines = data->nmachines;
    for (unsigned i = 0; i < nmachines; ++i) {
        if(sol->part[i].used == 0) {
            continue;
        }
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

int local_search_create_processing_list(solution *sol, local_search_data *data,
                                        int l) {
    int val = 0;
    Job *j1, *j2;
    int C;

    for (unsigned i = 0; i < data->nmachines; ++i) {
        C = 0;
        GPtrArray *machine = sol->part[i].machine;
        if(sol->part[i].used == 0) {
            continue;
        }

        for (unsigned j = 0; j < l; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            C += j1->processingime;
        }

        data->processing_list_f[i] = CC_SAFE_REALLOC(data->processing_list_f[i], data->njobs[i] - l, processing_list_data);


        for (unsigned j = 0; j < data->njobs[i] - l ; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            j2 = (Job *) g_ptr_array_index(machine, j + l);
            data->processing_list_f[i][j].pos = j;
            data->processing_list_f[i][j].p = C;
            C = C -  j1->processingime + j2->processingime;
        }

        qsort(data->processing_list_f[i], data->njobs[i] - l,
              sizeof(processing_list_data), compare_process_list);
    }

    return val;
}

int local_search_create_processing_list_b(solution *sol, local_search_data *data,
                                        int l) {
    int val = 0;
    Job *j1, *j2;
    int C;
    int njobs;

    for (unsigned i = 0; i < data->nmachines; ++i) {
        C = 0;
        njobs = data->njobs[i];
        GPtrArray *machine = sol->part[i].machine;
        if(sol->part[i].used == 0) {
            continue;
        }


        for (unsigned j = njobs - l; j < njobs; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            C += j1->processingime;
        }

        data->processing_list_b[i] = CC_SAFE_REALLOC(data->processing_list_b[i], njobs - l, processing_list_data);


        for (unsigned j = njobs - l; j > 0 ; --j) {
            j1 = (Job *) g_ptr_array_index(machine, j  - 1);
            j2 = (Job *) g_ptr_array_index(sol->part[i].machine, j + l - 1);
            data->processing_list_b[i][njobs - l - j].pos = j ;
            data->processing_list_b[i][njobs - l - j].p = C;
            C = C + j1->processingime - j2->processingime;
        }

        qsort(data->processing_list_b[i], njobs - l,
              sizeof(processing_list_data), compare_process_list_b);
    }

    return val;
}

int local_search_create_processing_list_swap(solution *sol, local_search_data *data,
                                        int l1, int l2) {
    int val = 0;
    Job *j1, *j2;
    int C;
    int njobs;
    int j;

    for (unsigned i = 0; i < data->nmachines; ++i) {
        njobs = data->njobs[i];
        GPtrArray *machine = sol->part[i].machine;
        if(sol->part[i].used == 0) {
            continue;
        }

        C = 0;
        for (unsigned j = l1; j < l1 + l2; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            C += j1->processingime;
        }

        data->processing_list_b[i] = CC_SAFE_REALLOC(data->processing_list_b[i], njobs - l1 - l2 +  1  , processing_list_data);

        for ( j = l1; j < njobs - l2  ; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            j2 = (Job *) g_ptr_array_index(machine, j + l2);
            data->processing_list_b[i][j - l1].pos = j;
            data->processing_list_b[i][j - l1].p = C;
            C = C - j1->processingime + j2->processingime;
        }
        data->processing_list_b[i][j - l1].pos = j;
        data->processing_list_b[i][j - l1].p = C;

        qsort(data->processing_list_b[i], njobs - l1 - l2 + 1 ,
              sizeof(processing_list_data), compare_process_list_b);

        C = 0;
        for (unsigned j = 0; j < l1 ; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            C += j1->processingime;
        }

        data->processing_list_f[i] = CC_SAFE_REALLOC(data->processing_list_f[i], njobs - l1 - l2 + 1, processing_list_data);

        for (unsigned j = 0; j < njobs - l1 - l2 + 1 ; ++j) {
            j1 = (Job *) g_ptr_array_index(machine, j);
            j2 = (Job *) g_ptr_array_index(machine, j + l1);
            data->processing_list_f[i][j].pos = j;
            data->processing_list_f[i][j].p = C;
            C = C - j1->processingime + j2->processingime;
        }

        qsort(data->processing_list_f[i], njobs - l1 - l2 + 1,
              sizeof(processing_list_data), compare_process_list);
    }

    return val;
}

int local_search_create_g(solution *sol, local_search_data *data) {
    int val = 0;
    int nmachines = sol->nmachines;
    Job *tmp;
    int t1, t2;
    int tw;
    int w;

    for (unsigned i = 0; i < nmachines; ++i) {
        if(sol->part[i].used == 0) {
            continue;
        } else {
            for(unsigned j = 0; j < sol->part[i].machine->len; ++j) {
                g_list_free_full(data->g[i][j], destroy_slope_t);
                data->g[i][j] = (GList *) NULL;
            }
            sol->part[i].used = 0;
        }


        int n_k = sol->part[i].machine->len;
        int P = 0;

        for (unsigned j = 0; j < n_k ; ++j) {
            GPtrArray *lateness_sort =  g_ptr_array_new();

            for (unsigned k = j; k < n_k; ++k) {
                g_ptr_array_add(lateness_sort, g_ptr_array_index(sol->part[i].machine, k));
            }

            g_ptr_array_sort_with_data(lateness_sort, local_search_compare_lateness,
                                       sol->c);
            tw = 0;
            w = 0;
            t1 = 0;
            int k;
            int move;

            move = 1;
            tmp = (Job *) g_ptr_array_index(lateness_sort, 0);
            for (k = 0  ; k < lateness_sort->len && move;) {
                move = tmp->weight * (sol->c[tmp->job] - P - tmp->duetime) > 0;

                if (move) {
                    tw += tmp->weight * (sol->c[tmp->job] - P - tmp->duetime);
                    w += tmp->weight;
                    k++;
                    tmp = (Job *) g_ptr_array_index(lateness_sort, k);
                }
            }


            t2 = tmp->duetime - sol->c[tmp->job] + P;
            local_search_add_slope_t(data, t1, t2, tw, w, i, j);


            for (unsigned l = k ; l < lateness_sort->len;) {
                tw = tw + w * (t2 - t1);
                t1 = t2;

                move = 1;
                tmp = (Job *) g_ptr_array_index(lateness_sort, l);

                while (move) {
                    w += tmp->weight;
                    l++;

                    if (l == lateness_sort->len) {
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

            g_ptr_array_free(lateness_sort, TRUE);
        }
    }

    return val;
}

void local_search_update_insertion(solution *sol, int i_best,
        int j_best, int k_best, int l) {
    Job *tmp;
    int old = sol->tw;
    sol->tw -= sol->part[k_best].tw;
    sol->part[k_best].c = 0;
    sol->part[k_best].tw = 0;

    for (unsigned i = 0; i < l; ++i) {
        tmp = (Job *) g_ptr_array_index(sol->part[k_best].machine, i_best);
        g_ptr_array_remove_index(sol->part[k_best].machine, i_best);
        g_ptr_array_insert(sol->part[k_best].machine, j_best, tmp);
    }

    for (unsigned i = 0; i < sol->part[k_best].machine->len; ++i) {
        tmp = (Job *) g_ptr_array_index(sol->part[k_best].machine, i);
        tmp->index = i;
        sol->part[k_best].c += tmp->processingime;
        sol->c[tmp->job] = sol->part[k_best].c;
        sol->part[k_best].tw += tmp->weight*CC_MAX(0, sol->c[tmp->job] - tmp->duetime);
    }


    sol->tw += sol->part[k_best].tw;
    printf("old - new:%d\n", old - sol->tw);
    sol->part[k_best].used = 1;
}

void local_search_update_swap(solution *sol, int i_best, int j_best, int k_best,
                              int l1, int l2) {
    Job *tmp;
    gpointer swap;
    int old = sol->tw;
    partlist *part = sol->part + k_best;
    sol->tw -= part->tw;
    part->c = 0;
    part->tw = 0;


    if(l1 == l2) {
        for(unsigned i = 0; i < l1; ++i) {
            CC_SWAP(g_ptr_array_index(part->machine,i_best + i), g_ptr_array_index(part->machine, j_best + i), swap);
        }
    } else if (l1 < l2){
        for(unsigned i = 0; i < l1; ++i) {
            CC_SWAP(g_ptr_array_index(part->machine,i_best + i), g_ptr_array_index(part->machine, j_best + i), swap);
        }

        for(unsigned i = l1; i < l2; ++i) {
            tmp = (Job*) g_ptr_array_index(part->machine, j_best + i);
            g_ptr_array_remove_index(part->machine, j_best + i);
            g_ptr_array_insert(part->machine, i_best + i, tmp);
        }
    } else {
        for(unsigned i = 0; i < l2; ++i) {
            CC_SWAP(g_ptr_array_index(part->machine,i_best + i), g_ptr_array_index(part->machine, j_best + i), swap);
        }

        for(unsigned i = l2; i < l1; ++i) {
            tmp = (Job*) g_ptr_array_index(part->machine, i_best + i);
            g_ptr_array_remove_index(part->machine, i_best + i);
            g_ptr_array_insert(part->machine, j_best + i - 1, tmp);
        }
    }

    for (unsigned i = 0; i < part->machine->len; ++i) {
        tmp = (Job *) g_ptr_array_index(sol->part[k_best].machine, i);
        tmp->index = i;
        part->c += tmp->processingime;
        sol->c[tmp->job] = part->c;
        part->tw += tmp->weight * CC_MAX(0, sol->c[tmp->job] - tmp->duetime);
    }

    sol->tw += part->tw;
    printf("old - new:%d\n", old - sol->tw);
    part->used = 1;
}

void local_search_forward_insertion(solution *sol, local_search_data *data,
                                    int l) {
    int pos, p;
    int **g, **h, *gg, **hh;
    int update;
    GList **iterators;
    Job *tmp_j;
    int max;
    int i_best, j_best, k_best;
    double runningtime = CCutil_zeit();

    g = CC_SAFE_MALLOC(sol->njobs, int *);
    h = CC_SAFE_MALLOC(sol->njobs, int *);
    hh = CC_SAFE_MALLOC(sol->njobs, int *);
    gg = CC_SAFE_MALLOC(sol->njobs, int);
    iterators = CC_SAFE_MALLOC(sol->njobs, GList *);

    for (unsigned i = 0; i < sol->njobs; ++i) {
        g[i] = CC_SAFE_MALLOC(sol->njobs, int);
        h[i] = CC_SAFE_MALLOC(sol->njobs, int);
        hh[i] = CC_SAFE_MALLOC(sol->njobs, int);
    }

    update = 0;
    max = 0;

    local_search_create_W(sol, data);
    local_search_create_g(sol, data);
    local_search_create_processing_list(sol, data, l);

    for (unsigned k = 0; k < sol->nmachines; ++k) {
        /** compute g */
        for (unsigned i = 0; i < data->njobs[k] - l; ++i) {
            pos = data->processing_list_f[k][i].pos;
            p = data->processing_list_f[k][i].p;
            GList *it = data->g[k][pos];

            for (unsigned j = pos + l; j < data->njobs[k]; ++j) {
                tmp_j = (Job *) g_ptr_array_index(sol->part[k].machine, j);
                int c = sol->c[tmp_j->job] - p;
                compute_it(&it, c);
                g[pos][j] = compute_g(&it, c);

            }
        }

        /** compute h */
        for (unsigned i = 0; i < data->njobs[k] - l; ++i) {
            GList *it = data->g[k][i + l];

            for (unsigned j = i + l; j < data->njobs[k]; ++j) {
                tmp_j = (Job *) g_ptr_array_index(sol->part[k].machine, j);
                int c = sol->c[tmp_j->job];
                compute_it(&it, c);
                h[i][j] = compute_g(&it, sol->c[tmp_j->job]);

            }
        }

        /** compute gg */
        for (unsigned i = 0; i < data->njobs[k] - l; ++i) {
            GList *it = data->g[k][i + l];
            int c;

            if (i != 0) {
                tmp_j = (Job *) g_ptr_array_index(sol->part[k].machine, i - 1);
                c  = sol->c[tmp_j->job];
            } else {
                c = 0;
            }

            compute_it(&it, c);
            gg[i] = compute_g(&it, c);
        }

        /** compute hh */
        for (unsigned j = l; j < data->njobs[k] - 1; ++j) {
            iterators[j] = data->g[k][j + 1];
        }

        iterators[data->njobs[k] - 1] = (GList *) NULL;

        for (unsigned i = 0; i < data->njobs[k] - l; ++i) {
            pos = data->processing_list_f[k][i].pos;
            p = data->processing_list_f[k][i].p;

            for (unsigned j = pos + l; j < data->njobs[k]; ++j) {
                if (iterators[j] != (GList *) NULL) {
                    tmp_j = (Job *) g_ptr_array_index(sol->part[k].machine, j);
                    int c = sol->c[tmp_j->job] - p ;
                    compute_it(&iterators[j], c);
                    hh[pos][j] = compute_g(&iterators[j], c);
                } else {
                    hh[pos][j] = 0;
                }

            }
        }

        for (unsigned i = 0; i < data->njobs[k] - l; ++i) {
            for (unsigned j = i + l; j < data->njobs[k]; ++j) {
                int tmp = 0;

                if (i != 0) {
                    tmp += data->W[k][i - 1];
                }

                tmp += g[i][j] - h[i][j];
                tmp += gg[i] - hh[i][j];
                tmp += data->W[k][data->njobs[k] - 1]  - data->W[k][j];

                if (sol->part[k].tw - tmp > max ) {
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
        solution_print(sol);
        local_search_update_insertion(sol, i_best, j_best, k_best, l);
        data->updated = 1;
        solution_print(sol);
    } else {
        data->updated = 0;
    }

    if(dbg_lvl() > 0) {
        printf("forward insertion with l = %d, running time = %f and improvement %d\n", l, CCutil_zeit() - runningtime, max);
    }

    for (unsigned i = 0; i < sol->njobs; ++i) {
        CC_IFFREE(g[i], int);
        CC_IFFREE(h[i], int);
        CC_IFFREE(hh[i], int);

    }

    CC_IFFREE(g, int *);
    CC_IFFREE(h, int *);
    CC_IFFREE(gg, int);
    CC_IFFREE(hh, int *);
    CC_IFFREE(iterators, GList *);
}

void local_search_backward_insertion(solution *sol, local_search_data *data,
                                    int l) {
    int p;
    int c;
    int **g, **h, *gg, **hh;
    int update;
    GList **iterators;
    Job *tmp_j;
    int max;
    int i_best, j_best, k_best;
    double runningtime = CCutil_zeit();

    g = CC_SAFE_MALLOC(sol->njobs, int *);
    h = CC_SAFE_MALLOC(sol->njobs, int *);
    hh = CC_SAFE_MALLOC(sol->njobs, int *);
    iterators = CC_SAFE_MALLOC(sol->njobs, GList *);
    gg = CC_SAFE_MALLOC(sol->njobs, int);

    for (unsigned i = 0; i < sol->njobs; ++i) {
        g[i] = CC_SAFE_MALLOC(sol->njobs, int);
        h[i] = CC_SAFE_MALLOC(sol->njobs, int);
        hh[i] = CC_SAFE_MALLOC(sol->njobs, int);
    }

    update = 0;
    max = 0;

    local_search_create_processing_list_b(sol, data, 2);
    local_search_create_W(sol, data);
    local_search_create_g(sol, data);

    for (unsigned k = 0; k < sol->nmachines; ++k) {
        int njobs = data->njobs[k];
        GPtrArray *machine = sol->part[k].machine;
        /** compute g */
        for (int i = njobs - l; i > 0; --i) {
            GList *it = data->g[k][i];

            for (int j = 0; j < i ; ++j) {
                if(j == 0) {
                    c = 0;
                } else {
                    tmp_j = (Job *) g_ptr_array_index(machine, j - 1);
                    c = sol->c[tmp_j->job];
                }
                compute_it(&it, c);
                g[i][j] = compute_g(&it, c);

            }
        }

        /** compute h */
        for(unsigned j = 0; j < njobs - l; ++j) {
            h[njobs - l][j] = 0;
        }

        p = 0;
        for(unsigned i = njobs - l - 1; i < njobs - 1; ++i) {
            tmp_j = (Job *) g_ptr_array_index(machine, i);
            p += tmp_j->processingime;
        }

        for (int i = njobs - l - 1; i > 0; --i) {
            GList *it = data->g[k][i + l];

            for (int j = 0 ; j < i ; ++j) {
                if(j == 0) {
                    c = p;
                } else {
                    tmp_j = (Job *) g_ptr_array_index(machine, j - 1);
                    c = p + sol->c[tmp_j->job];
                }
                compute_it(&it, c);
                h[i][j] = compute_g(&it, c );
            }

            tmp_j = (Job *) g_ptr_array_index(machine, i + l - 1);
            p -= tmp_j->processingime;
            tmp_j = (Job *) g_ptr_array_index(machine, i - 1);
            p += tmp_j->processingime;
        }

        /** compute gg */
        for (unsigned i = njobs - l; i > 0; --i) {
            GList *it = data->g[k][i];
            tmp_j = (Job *) g_ptr_array_index(machine, i + l - 1);
            c  = sol->c[tmp_j->job];
            compute_it(&it, c);
            gg[i] = compute_g(&it, c);
        }

        /** compute hh */
        for (unsigned j = 0; j < njobs - l; ++j) {
            iterators[j] = data->g[k][j];
        }

        for (unsigned i = njobs - l - 1; i > 0;--i) {
            p = data->processing_list_b[k][njobs -l - i].p;
            int pos = data->processing_list_b[k][njobs - l - i].pos;

            for (unsigned j = 0; j < pos; ++j) {
                    if(j == 0) {
                        c = p;
                    } else {
                        tmp_j = (Job *) g_ptr_array_index(machine, j - 1);
                        c = p + sol->c[tmp_j->job];
                    }
                    compute_it(&iterators[j], c);
                    hh[pos][j] = compute_g(&iterators[j], c);
            }
        }

        for (unsigned i = njobs - l; i > 0 ; --i) {
            for (int j = i - 1; j >= 0; --j) {
                int t = 0;

                if (j != 0) {
                    t += data->W[k][j - 1];
                }

                t += g[i][j] - h[i][j];
                t += hh[i][j] - gg[i];
                t += data->W[k][njobs - 1]  - data->W[k][i + l - 1];

                if (sol->part[k].tw - t > max ) {
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
        solution_print(sol);
        local_search_update_insertion(sol, i_best, j_best, k_best, l);
        solution_print(sol);
    }

    if(dbg_lvl() > 0) {
        printf("backward insertion with l = %d, running time = %f and improvement %d\n", l, CCutil_zeit() - runningtime, max);
    }

    for (unsigned i = 0; i < sol->njobs; ++i) {
        CC_IFFREE(g[i], int);
        CC_IFFREE(h[i], int);
        CC_IFFREE(hh[i], int);

    }

    CC_IFFREE(g, int *);
    CC_IFFREE(h, int *);
    CC_IFFREE(gg, int);
    CC_IFFREE(hh, int *);
    CC_IFFREE(iterators, GList *);
}

void local_search_swap_intra(solution *sol, local_search_data *data,
                                    int l1, int l2) {
    int pos, p;
    int **B2_1, **B2_2, **B3_1, **B3_2, **B4_1, **B4_2;
    int update;
    GList **iterators;
    Job *tmp_j;
    int max;
    int i_best = -1, j_best = -1, k_best = -1;
    double runningtime = CCutil_zeit();

    B2_1 = CC_SAFE_MALLOC(sol->njobs, int *);
    B2_2 = CC_SAFE_MALLOC(sol->njobs, int *);
    B3_1 = CC_SAFE_MALLOC(sol->njobs, int*);
    B3_2 = CC_SAFE_MALLOC(sol->njobs, int *);
    B4_1 = CC_SAFE_MALLOC(sol->njobs, int *);
    B4_2 = CC_SAFE_MALLOC(sol->njobs, int *);
    iterators = CC_SAFE_MALLOC(sol->njobs, GList *);


    for (unsigned i = 0; i < sol->njobs; ++i) {
        B2_1[i] = CC_SAFE_MALLOC(sol->njobs, int);
        B2_2[i] = CC_SAFE_MALLOC(sol->njobs, int);
        B3_1[i] = CC_SAFE_MALLOC(sol->njobs, int);
        B3_2[i] = CC_SAFE_MALLOC(sol->njobs, int);
        B4_1[i] = CC_SAFE_MALLOC(sol->njobs, int);
        B4_2[i] = CC_SAFE_MALLOC(sol->njobs, int);
    }

    update = 0;
    max = 0;

    local_search_create_W(sol, data);
    local_search_create_processing_list_swap(sol, data, l1, l2);
    local_search_create_g(sol, data);

    for (unsigned k = 0; k < sol->nmachines; ++k) {
        int njobs = sol->part[k].machine->len;
        GPtrArray *machine = sol->part[k].machine;
        /** compute g */
        for (unsigned i = 0; i < njobs - l1 - l2 + 1; ++i) {
            pos = data->processing_list_f[k][i].pos;
            p = data->processing_list_f[k][i].p;
            GList *it = data->g[k][pos];

            for (unsigned j = pos + l1; j < njobs - l2 + 1 ; ++j) {
                tmp_j = (Job *) g_ptr_array_index(machine, j + l2 - 1);
                int c = sol->c[tmp_j->job] - p; 
                compute_it(&it, c);
                B2_1[pos][j] = compute_g(&it, c);
            }
        }

        /** compute h */
        for (unsigned i = 0; i < njobs - l1 - l2 + 1; ++i) {
            GList *it = data->g[k][i + l1];

            for (unsigned j = i + l1; j < njobs - l2 + 1 ; ++j) {
                tmp_j = (Job *) g_ptr_array_index(machine, j + l2 - 1);
                int c = sol->c[tmp_j->job]; 
                compute_it(&it, c);
                B2_2[i][j] = compute_g(&it, c);
            }
        }


        /** compute B3_1 */
        for(unsigned i = 0; i < njobs - l1 - l2 + 1 ; ++i) {
            iterators[i] = data->g[k][i + l1];
        }

        for(unsigned j = l1; j < njobs - l2 + 1 ; ++j) {
            pos = data->processing_list_b[k][j - l1].pos;
            p = data->processing_list_b[k][j - l1].p;
            for(unsigned i = 0; i < pos - l1 + 1 ; ++i) {
                int c = p;
                if(i != 0) {
                    tmp_j = (Job *) g_ptr_array_index(machine, i - 1);
                    c += sol->c[tmp_j->job];
                }
                compute_it(&iterators[i], c);
                B3_1[i][pos] = compute_g(&iterators[i], c);
            }
        }

        /** compute B3_2 */
        for(unsigned j = l1; j < njobs - l2 + 1 ; ++j) {
            iterators[j] = data->g[k][j];
        }

        for(unsigned i = 0; i < njobs - l1 - l2 + 1; ++i) {
            pos = data->processing_list_f[k][i].pos;
            p = data->processing_list_f[k][i].p;
            for(unsigned j = pos + l1; j < njobs - l2 + 1 ; ++j) {
                tmp_j = (Job *) g_ptr_array_index(machine, j + l2 - 1);
                int c = sol->c[tmp_j->job] - p;
                compute_it(&iterators[j], c);
                B3_2[pos][j] = compute_g(&iterators[j], c);
            }
        }

        /** compute B4_1 */
        for(unsigned j = l1; j < njobs - l2 + 1 ; ++j) {
            GList *it = data->g[k][j];
            for(unsigned i = 0; i < j - l1 + 1; ++i) {
                int c = 0;
                if(i != 0) {
                    c = sol ->c[((Job *) g_ptr_array_index(machine, i - 1))->job];
                }
                compute_it(&it, c);
                B4_1[i][j] = compute_g(&it, c);
            }
        }

        /** compute B4_2 */
        for(unsigned j = l1; j < njobs - l2 + 1 ; ++j) {
            pos = data->processing_list_b[k][j - l1].pos;
            p = data->processing_list_b[k][j - l1].p;
            GList *it;
            if(pos + l2 == njobs) {
                it = NULL;
            } else {
                it = data->g[k][pos + l2];
            }
            for(unsigned i = 0; i < pos - l1 + 1; ++i) {
                int c = p;
                if(i != 0) {
                    c += sol->c[((Job *)g_ptr_array_index(machine, i - 1))->job];
                }
                compute_it(&it, c);
                if(it != NULL) {
                    B4_2[i][pos] = compute_g(&it,c);
                } else {
                    B4_2[i][pos] = 0;
                }
            }
        }

        for (unsigned i = 0; i < njobs - l1 - l2 + 1  ; ++i) {
            for (unsigned j = i + l1; j < njobs - l2 + 1 ; ++j) {
                int t = 0;

                if (i != 0) {
                    t += data->W[k][i - 1];
                }

                t += B2_1[i][j] - B2_2[i][j];
                t += B3_1[i][j] - B3_2[i][j];
                t += B4_1[i][j] - B4_2[i][j];
                t += data->W[k][njobs - 1]  - data->W[k][j + l2 - 1];

                if (sol->part[k].tw - t > max ) {
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
    solution_print(sol);
    if (update) {
        local_search_update_swap(sol, i_best, j_best, k_best, l1, l2);
        data->updated = 1;
        solution_print(sol);
    } else {
        data->updated = 0;
    }

    if(dbg_lvl() > 0) {
        printf("intra swap with l1 = %d and l2 = %d, running time = %f and improvement %d on machine %d on places %d %d\n", l1,l2, CCutil_zeit() - runningtime, max, k_best,i_best,j_best);
    }

    for (unsigned i = 0; i < sol->njobs; ++i) {
        CC_IFFREE(B2_1[i], int);
        CC_IFFREE(B2_2[i], int);
        CC_IFFREE(B3_1[i], int);
        CC_IFFREE(B3_2[i], int);
        CC_IFFREE(B4_1[i], int);
        CC_IFFREE(B4_2[i], int);

    }

    CC_IFFREE(B2_1, int *);
    CC_IFFREE(B2_2, int *);
    CC_IFFREE(B3_1, int*);
    CC_IFFREE(B3_2, int *);
    CC_IFFREE(B4_1, int *);
    CC_IFFREE(B4_2, int *);
    CC_IFFREE(iterators, GList *);
}