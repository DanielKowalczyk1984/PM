#include <interval.h>
#include <job.h>
#include <solution.h>
#include <stdlib.h>
#include <string.h>
#include <util.h>

gint g_sort_jobs_key(const void* a, const void* b, void* data);
gint g_sort_jobs_weight(gconstpointer a, gconstpointer b, void* data);

void solution_init(Solution* sol) {
    if (sol) {
        sol->part = (PartList*)NULL;
        sol->perm = (Job**)NULL;
        sol->c = (int*)NULL;
        sol->u = (int*)NULL;
        sol->nb_machines = 0;
        sol->nb_jobs = 0;
        sol->nb_intervals = 0;
        sol->tw = 0;
        sol->b = 0;
        sol->off = 0;
    }
}

void solution_free(Solution** sol) {
    if (*sol) {
        for (int i = 0; i < (*sol)->nb_machines; ++i) {
            for(int j = 0; j < (*sol)->nb_intervals; ++j) {
                g_ptr_array_free((*sol)->part[i].Q[j],TRUE);
                g_ptr_array_free((*sol)->part[i].Q_in[j],TRUE);
            }
            partlist_free((*sol)->part + i);
        }

        CC_IFFREE((*sol)->part, PartList);
        CC_IFFREE((*sol)->perm, Job*);
        CC_IFFREE((*sol)->c, int);
        CC_IFFREE((*sol)->u, int);
        CC_IFFREE((*sol)->u_in, int);
        CC_IFFREE((*sol), Solution);
    }
}

Solution* solution_alloc(int nb_interval,int nb_machines, int nb_jobs, int off) {
    int       val = 0;
    int       i;
    Solution* sol = CC_SAFE_MALLOC(1, Solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory");
    solution_init(sol);
    sol->nb_machines = nb_machines;
    sol->nb_jobs = nb_jobs;
    sol->nb_intervals = nb_interval;
    sol->tw = 0;
    sol->b = 0;
    sol->off = off;
    sol->part = CC_SAFE_MALLOC(nb_machines, PartList);
    CCcheck_NULL_2(sol->part, "Failed to allocate memory to part");

    for (i = 0; i < nb_machines; ++i) {
        partlist_init(sol->part + i);
        sol->part[i].Q = CC_SAFE_MALLOC(nb_interval, GPtrArray*);
        sol->part[i].Q_in = CC_SAFE_MALLOC(nb_interval, GPtrArray*);
        for(int j = 0; j < nb_interval; j++) {
            sol->part[i].Q[j] = g_ptr_array_new();
            sol->part[i].Q_in[j] = g_ptr_array_new();
        }
    }

    sol->perm = CC_SAFE_MALLOC(nb_jobs, Job*);
    CCcheck_NULL_2(sol->perm, "Failed to allocate memory to perm");
    sol->c = CC_SAFE_MALLOC(nb_jobs, int);
    CCcheck_NULL_2(sol->c, "Failed to allocate memory");
    fill_int(sol->c, sol->nb_jobs, 0);
    sol->u = CC_SAFE_MALLOC(nb_jobs, int);
    CCcheck_NULL_2(sol->u, "Failed to allocate memory")
    fill_int(sol->u, nb_jobs, 0);
    sol->u_in = CC_SAFE_MALLOC(nb_jobs, int);
    CCcheck_NULL_2(sol->u_in, "Failed to allocate memory");
    fill_int(sol->u_in, nb_jobs, -1);

    for (i = 0; i < nb_jobs; ++i) {
        sol->perm[i] = (Job*)NULL;
    }

CLEAN:

    if (val) {
        solution_free(&sol);
    }

    return sol;
}

gint g_sort_jobs_key(const void* a, const void* b, void* data) {
    (void)data;
    const int* v = &(((const Job*)a)->job);
    const int* w = &(((const Job*)b)->job);
    return *v - *w;
}

gint g_sort_jobs_weight(gconstpointer a, gconstpointer b, void* data) {
    (void)data;
    const int* v = &(((const Job*)a)->weight);
    const int* w = &(((const Job*)b)->weight);
    return -(*v - *w);
}

static void g_print_jobs(gpointer j, gpointer data) {
    Job* tmp = (Job*)j;
    printf("%d ", tmp->job);
}

void solution_print(Solution* sol) {
    for (int i = 0; i < sol->nb_machines; ++i) {
        printf("Machine %-1d: ", i);
        g_ptr_array_foreach(sol->part[i].machine, g_print_jobs, sol);
        printf("with C =  %d, wC = %d and %u jobs\n", sol->part[i].c,
               sol->part[i].tw, sol->part[i].machine->len);
    }

    printf("with total weighted tardiness %d\n", sol->tw + sol->off);
}

int solution_copy(Solution* dest, Solution* src) {
    int val = 0;
    dest = solution_alloc(src->nb_intervals, src->nb_machines, src->nb_jobs, src->off);
    CCcheck_val_2(val, "Failed in  solution_alloc");
    dest->tw = src->tw;
    dest->b = src->b;
    dest->off = src->off;

    for (int i = 0; i < dest->nb_machines; i++) {
        dest->part[i].tw = src->part[i].tw;
        dest->part[i].c = src->part[i].c;
    }

CLEAN:

    if (val) {
        solution_free(&dest);
    }

    return val;
}

int solution_update(Solution* dest, Solution* src) {
    int val = 0;
    dest->tw = src->tw;
    dest->b = src->b;
    dest->nb_machines = src->nb_machines;
    dest->nb_jobs = src->nb_jobs;
    dest->off = src->off;

    for (int i = 0; i < dest->nb_machines; i++) {
        g_ptr_array_remove_range(dest->part[i].machine, 0,
                                 dest->part[i].machine->len);

        for (unsigned j = 0; j < src->part[i].machine->len; ++j) {
            g_ptr_array_add(dest->part[i].machine,
                            g_ptr_array_index(src->part[i].machine, j));
        }

        dest->part[i].tw = src->part[i].tw;
        dest->part[i].c = src->part[i].c;

        for(int j = 0; j < src->nb_intervals; j++) {
            g_ptr_array_remove_range(dest->part[i].Q[j], 0, dest->part[i].Q[j]->len);
            g_ptr_array_remove_range(dest->part[i].Q_in[j], 0, dest->part[i].Q_in[j]->len);
            for(int k = 0; k < src->part[i].Q[j]->len;k++) {
                g_ptr_array_add(dest->part[i].Q[j], g_ptr_array_index(src->part[i].Q[j],k));
            }
            for(int k = 0; k < src->part[i].Q_in[j]->len;k++) {
                g_ptr_array_add(dest->part[i].Q_in[j], g_ptr_array_index(src->part[i].Q_in[j],k));
            }
        }
    }

    memcpy(dest->perm, src->perm, src->nb_jobs * sizeof(Job*));
    memcpy(dest->c, src->c, dest->nb_jobs * sizeof(int));
    memcpy(dest->u, src->c, dest->nb_jobs * sizeof(int));
    return val;
}

void g_reset_num_layers(gpointer data, gpointer user_data) {
    Job* j = (Job*)data;
    j->num_layers = 0;
}

void reset_nb_layers(GPtrArray* jobs) {
    g_ptr_array_foreach(jobs, g_reset_num_layers, NULL);
}

void g_set_sol_perm(gpointer data, gpointer user_data) {
    Job*      j = (Job*)data;
    Solution* sol = (Solution*)user_data;
    sol->perm[j->job] = j;
}

void solution_calculate_machine(Solution* sol, int m) {
    if (m < sol->nb_machines) {
        PartList*  part = sol->part + m;
        GPtrArray* machine = sol->part[m].machine;
        sol->tw -= part->tw;
        part->tw = 0;
        part->c = 0;

        for (unsigned i = 0; i < machine->len; ++i) {
            Job* tmp = (Job*)g_ptr_array_index(machine, i);
            part->c += tmp->processing_time;
            sol->c[tmp->job] = part->c;
            part->tw += value_Fj(sol->c[tmp->job], tmp);
        }

        sol->tw += part->tw;
    }
}

void solution_calculate_all(Solution* sol) {
    for (int i = 0; i < sol->nb_machines; ++i) {
        solution_calculate_machine(sol, i);
    }
}

void solution_calculate_partition_machine(Solution* sol, GPtrArray* intervals,
                                          int m) {
    if (m < sol->nb_machines) {
        GPtrArray* machine = sol->part[m].machine;
        for(int i = 0; i < sol->nb_intervals;i++) {
            g_ptr_array_free(sol->part[m].Q[i],TRUE);
            sol->part[m].Q[i] = g_ptr_array_new();
            g_ptr_array_free(sol->part[m].Q_in[i],TRUE);
            sol->part[m].Q_in[i] = g_ptr_array_new();
        }
        int        iter = 0;

        for (unsigned i = 0; i < machine->len; ++i) {
            Job*      tmp = (Job*)g_ptr_array_index(machine, i);
            interval* I = (interval*)g_ptr_array_index(intervals, iter);
            while (!(sol->c[tmp->job] <= I->b && I->a < sol->c[tmp->job])) {
                iter++;
                I = (interval*)g_ptr_array_index(intervals, iter);
            }
            sol->u[tmp->job] = iter; 
            g_ptr_array_add(sol->part[m].Q[iter], tmp);
            if(sol->c[tmp->job] - tmp->processing_time > I->a || (iter == 0 && sol->c[tmp->job] - tmp->processing_time == 0 )) {
                g_ptr_array_add(sol->part[m].Q_in[iter], tmp);
                sol->u_in[tmp->job] = iter;
            } 
        }
    }
}

void solution_calculate_partition_all(Solution* sol, GPtrArray* intervals) {
    for (int i = 0; i < sol->nb_machines; ++i) {
        solution_calculate_partition_machine(sol, intervals, i);
    }
}

// static void calculate_partition(Solution* sol, GPtrArray* intervals, int m,
//                                 int* u, int* last) {
//     int   count = 0;
//     int   cur = *last;
//     void* tmp;

//     GPtrArray* machine = sol->part[m].machine;
//     interval*  I = (interval*)g_ptr_array_index(intervals, *u);
//     Job*       first_job = g_ptr_array_index(machine, 0);
//     Job*       i = (Job*)g_ptr_array_index(machine, cur);
//     Job*       j;

//     while (sol->c[i->job] > I->a && sol->c[i->job] <= I->b &&
//            sol->c[i->job] - i->processing_time <= I->b &&
//            sol->c[i->job] - i->processing_time > I->a) {
//         if (first_job != i) {
//             count++;
//             i = (Job*)g_ptr_array_index(machine, --cur);
//         } else {
//             count++;
//             break;
//         }
//     }

//     if (count > 0) {
//         if (first_job->job != i->job) {
//             g_qsort_with_data(machine->pdata + cur + 1, count, sizeof(Job*),
//                               compare_interval, I);
//             j = (Job*)g_ptr_array_index(machine, cur + 1);
//         } else {
//             g_qsort_with_data(machine->pdata, count + 1, sizeof(Job*),
//                               compare_interval, I);
//             i = (Job*)g_ptr_array_index(machine, 0);
//             j = (Job*)g_ptr_array_index(machine, 1);
//         }

//         if (sol->c[i->job] > I->a && sol->c[i->job] <= I->b) {
//             if (compare_interval(&i, &j, I) < 0) {
//                 cur--;
//                 *last = cur;
//                 if (cur >= 0) {
//                     i = (Job*)g_ptr_array_index(machine, cur);
//                     *u = CC_MIN(*u - 1, sol->u[i->job]);
//                 } else {
//                     *u = -1;
//                 }
//             } else {
//                 sol->c[j->job] =
//                     sol->c[i->job] - i->processing_time + j->processing_time;
//                 sol->c[i->job] = sol->c[j->job] + i->processing_time;
//                 CC_SWAP(g_ptr_array_index(machine, cur),
//                         g_ptr_array_index(machine, cur + 1), tmp);
//                 if (sol->c[j->job] > I->a && sol->c[j->job] <= I->b) {
//                     g_qsort_with_data(machine->pdata + cur, count + 1,
//                                       sizeof(Job*), compare_interval, I);
//                     *last = cur;
//                     if (cur >= 0) {
//                         i = (Job*)g_ptr_array_index(machine, cur);
//                         *u = CC_MIN(*u - 1, sol->u[i->job]);
//                     } else {
//                         *u = -1;
//                     }
//                 }
//             }
//         }
//     } else {
//         cur--;
//         *last = cur;
//         if (cur >= 0) {
//             i = (Job*)g_ptr_array_index(machine, cur);
//             *u = CC_MIN(*u - 1, sol->u[i->job]);
//         } else {
//             *u = -1;
//         }
//     }

//     solution_calculate_machine(sol, m);
// }
static int next_interval_reversed(int u, PartList* part) {
    u--;
    while(u != -1) {
        GPtrArray *Q = part->Q_in[u];
        if(Q->len > 0) return u;
        u--;
    }

    return -1;

}
int solution_canonical_order(Solution* sol, GPtrArray* intervals) {
    int val = 0;

    solution_calculate_partition_all(sol, intervals);
    sol->tw = 0;

    for (int it = 0; it < sol->nb_machines; ++it) {
        PartList *part = sol->part + it;
        GPtrArray* machine = part->machine;
        int        last = machine->len - 1;
        Job*       i = (Job*)g_ptr_array_index(machine, last);
        int        u = sol->u[i->job];
        while (u != -1) {
            interval *I = (interval *) g_ptr_array_index(intervals, u);
            GPtrArray *Q = part->Q[u];
            GPtrArray *Q_in = part->Q_in[u];
            if(Q_in->len > 0) {
                if(Q_in->len + 1 == Q->len) {
                    Job *first = (Job *) g_ptr_array_index(Q_in, 0);
                    int C = sol->c[first->job] - first->processing_time;
                    #ifndef NDEBUG
                    Job *last = (Job *) g_ptr_array_index(Q_in, Q_in->len - 1);
                    int C_last = sol->c[last->job];
                    #endif
                    g_ptr_array_sort_with_data(Q_in, g_compare_interval_data, I);
                    for(int j = 0; j < Q_in->len;j++) {
                        Job *tmp = g_ptr_array_index(Q_in, j);
                        g_ptr_array_index(Q, j + 1) = tmp;
                        C += tmp->processing_time;
                        sol->c[tmp->job] = C;
                    }
                    assert(C == C_last);
                    Job *tmp_out = (Job*) g_ptr_array_index(Q, 0);
                    Job *tmp_in = (Job*) g_ptr_array_index(Q_in, 0);
                    if(g_compare_interval_data(&tmp_out, &tmp_in,I) < 0) {
                        assert(sol->c[tmp_in->job] - tmp_in->processing_time == sol->c[tmp_out->job]);
                        u = next_interval_reversed(u, part);
                    } else {
                        Job *tmp = NULL;
                        CC_SWAP(g_ptr_array_index(Q, 0), g_ptr_array_index(Q, 1), tmp);
                        g_ptr_array_index(Q_in, 0) = tmp_out;
                        sol->c[tmp_in->job] = sol->c[tmp_out->job] - tmp_out->processing_time + tmp_in->processing_time;
                        sol->c[tmp_out->job] = sol->c[tmp_in->job] + tmp_out->processing_time;

                        if(sol->c[tmp_in->job] <= I->b && sol->c[tmp_in->job] > I->a) {
                            sol->u[tmp_out->job] = u;
                            sol->u_in[tmp_out->job] = u;
                            sol->u[tmp_in->job] = u;
                            int C = sol->c[tmp_in->job];
                            assert(C == sol->c[tmp_out->job] - tmp_out->processing_time);
                            g_ptr_array_sort_with_data(Q_in, g_compare_interval_data, I);
                            for(int j = 0; j < Q_in->len;j++) {
                                Job *tmp = g_ptr_array_index(Q_in, j);
                                g_ptr_array_index(Q, j + 1) = tmp;
                                C += tmp->processing_time;
                                sol->c[tmp->job] = C;
                            }
                            u = next_interval_reversed(u, part);
                        } else {
                            if(sol->c[tmp_out->job] <= I->b && sol->c[tmp_out->job] - tmp_out->processing_time > I->a) {
                                sol->u[tmp_out->job] = u;
                                sol->u_in[tmp_out->job] = u; 
                            } else {
                                // assert(sol->c[tmp_out->job] <= I->b && sol->c[tmp_out->job] > I->a);
                                sol->u[tmp_out->job] = u;
                                sol->u[tmp_out->job] = -1;
                                g_ptr_array_remove(Q_in, tmp_out);
                            }
                            g_ptr_array_remove(Q, tmp_in);
                            int old_u = u - 1;
                            while(old_u != -1) {
                                interval* tmp_I = (interval*)g_ptr_array_index(intervals, old_u);
                                if(sol->c[tmp_in->job] <= tmp_I->b && sol->c[tmp_in->job] > tmp_I->a) {
                                    g_ptr_array_add(part->Q[old_u], tmp_in);
                                    sol->u[tmp_in->job] = old_u;
                                    if(sol->c[tmp_in->job] - tmp_in->processing_time > tmp_I->a) {
                                        g_ptr_array_add(part->Q_in[old_u], tmp_in);
                                        sol->u_in[tmp_in->job] = old_u;
                                    }
                                    break;
                                }
                                old_u--;
                            }
                        }
                    }
                } else {
                    assert(u == 0);
                    g_qsort_with_data(Q->pdata, Q->len, sizeof(Job*), g_compare_interval_data, I);
                    int C = 0;
                    for(int j = 0; j < Q->len;j++) {
                        Job *tmp = g_ptr_array_index(Q, j);
                        C += tmp->processing_time;
                        sol->c[tmp->job] = C;
                    }
                    u--;
                }
            } else {
                u--;
            }
            // calculate_partition(sol, intervals, it, &u, &last);
        }

        g_ptr_array_remove_range(machine, 0, machine->len);
        part->tw = 0;
        part->c = 0;

        for(int uu = 0; uu < intervals->len; uu++) {
            GPtrArray *Q = part->Q[uu];
            if(Q->len > 0) {
                for(int k = 0; k < Q->len; k++) {
                    Job *tmp = g_ptr_array_index(Q, k);
                    g_ptr_array_add(machine, g_ptr_array_index(Q, k));
                    part->c += tmp->processing_time;
                    assert(part->c == sol->c[tmp->job]);
                    sol->c[tmp->job] = part->c;
                    part->tw += value_Fj(part->c, tmp);
                }
            }
        }

        sol->tw += part->tw;
        // solution_calculate_machine(sol, it);
    }

    return val;
}

int solution_arctime_order(Solution* sol) {
    int val = 0;
    solution_print(sol);

    for (int it = 0; it < sol->nb_machines; it++) {
        GPtrArray* tmp = sol->part[it].machine;
        int        C = ((Job*)g_ptr_array_index(tmp, 0))->processing_time;
        int        i = 0;
        int        j = 1;

        while (i < tmp->len - 2 && j < tmp->len - 1) {
            Job* tmp_i = (Job*)g_ptr_array_index(tmp, i);
            Job* tmp_j = (Job*)g_ptr_array_index(tmp, j);
            if (tmp_i->job < tmp_j->job) {
                if (arctime_diff_Fij(C, tmp_i, tmp_j) >= 0) {
                    void* tmp_job = NULL;
                    CC_SWAP(g_ptr_array_index(tmp, i),
                            g_ptr_array_index(tmp, j), tmp_job);
                    C = ((Job*)g_ptr_array_index(tmp, 0))->processing_time;
                    i = 0;
                    j = 1;
                } else {
                    C += tmp_j->processing_time;
                    i++;
                    j++;
                }
            } else {
                if (arctime_diff_Fij(C, tmp_i, tmp_j) > 0) {
                    void* tmp_job = NULL;
                    CC_SWAP(g_ptr_array_index(tmp, i),
                            g_ptr_array_index(tmp, j), tmp_job);
                    C = ((Job*)g_ptr_array_index(tmp, 0))->processing_time;
                    i = 0;
                    j = 1;
                } else {
                    C += tmp_j->processing_time;
                    i++;
                    j++;
                }
            }
        }
    }
    solution_calculate_all(sol);

    solution_print(sol);
    return val;
}
