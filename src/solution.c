#include <interval.h>
#include <job.h>
#include <solution.h>
#include <stdlib.h>
#include <string.h>
#include <util.h>

gint comparefunc(const void* a, const void* b, void* data);
gint compare_func(const void* a, const void* b);
gint order_weight(gconstpointer a, gconstpointer b, void* data);

void solution_init(Solution* sol) {
    if (sol) {
        sol->part = (PartList*)NULL;
        sol->perm = (Job**)NULL;
        sol->c = (int*)NULL;
        sol->u = (int*)NULL;
        sol->nb_machines = 0;
        sol->nb_jobs = 0;
        sol->tw = 0;
        sol->b = 0;
        sol->off = 0;
    }
}

void solution_free(Solution** sol) {
    if (*sol) {
        for (int i = 0; i < (*sol)->nb_machines; ++i) {
            partlist_free((*sol)->part + i);
        }

        CC_IFFREE((*sol)->part, PartList);
        CC_IFFREE((*sol)->perm, Job*);
        CC_IFFREE((*sol)->c, int);
        CC_IFFREE((*sol)->u, int);
        CC_IFFREE((*sol), Solution);
    }
}

Solution* solution_alloc(int nb_machines, int nb_jobs, int off) {
    int       val = 0;
    int       i;
    Solution* sol = CC_SAFE_MALLOC(1, Solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory");
    solution_init(sol);
    sol->nb_machines = nb_machines;
    sol->nb_jobs = nb_jobs;
    sol->tw = 0;
    sol->b = 0;
    sol->off = off;
    sol->part = CC_SAFE_MALLOC(nb_machines, PartList);
    CCcheck_NULL_2(sol->part, "Failed to allocate memory to part");

    for (i = 0; i < nb_machines; ++i) {
        partlist_init(sol->part + i);
        (sol->part + i)->key = i;
    }

    sol->perm = CC_SAFE_MALLOC(nb_jobs, Job*);
    CCcheck_NULL_2(sol->perm, "Failed to allocate memory to perm");
    sol->c = CC_SAFE_MALLOC(nb_jobs, int);
    CCcheck_NULL_2(sol->c, "Failed to allocate memory");
    fill_int(sol->c, sol->nb_jobs, 0);
    sol->u = CC_SAFE_MALLOC(nb_jobs, int);
    CCcheck_NULL_2(sol->u, "Failed to allocate memory")
        fill_int(sol->u, nb_jobs, 0);

    for (i = 0; i < nb_jobs; ++i) {
        sol->perm[i] = (Job*)NULL;
    }

CLEAN:

    if (val) {
        solution_free(&sol);
    }

    return sol;
}

gint comparefunc(const void* a, const void* b, void* data) {
    (void)data;
    const int* v = &(((const Job*)a)->job);
    const int* w = &(((const Job*)b)->job);
    return *v - *w;
}

gint order_weight(gconstpointer a, gconstpointer b, void* data) {
    (void)data;
    const int* v = &(((const Job*)a)->weight);
    const int* w = &(((const Job*)b)->weight);
    return -(*v - *w);
}

static void print_machine(gpointer j, gpointer data) {
    Job* tmp = (Job*)j;
    int  C = ((Solution*)data)->c[tmp->job];
    printf("%d (%d) ", tmp->job, C);
}

void solution_print(Solution* sol) {
    for (int i = 0; i < sol->nb_machines; ++i) {
        printf("Machine %-1d: ", sol->part[i].key);
        g_ptr_array_foreach(sol->part[i].machine, print_machine, sol);
        printf("with C =  %d, wC = %d and %u jobs\n", sol->part[i].c,
               sol->part[i].tw, sol->part[i].machine->len);
    }

    printf("with total weighted tardiness %d\n", sol->tw + sol->off);
}

int solution_copy(Solution* dest, Solution* src) {
    int val = 0;
    dest = solution_alloc(src->nb_machines, src->nb_jobs, src->off);
    CCcheck_val_2(val, "Failed in  solution_alloc");
    dest->tw = src->tw;
    dest->b = src->b;
    dest->off = src->off;

    for (int i = 0; i < dest->nb_machines; i++) {
        dest->part[i].key = src->part[i].key;
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
    }

    memcpy(dest->perm, src->perm, src->nb_jobs * sizeof(Job*));
    memcpy(dest->c, src->c, dest->nb_jobs * sizeof(int));
    memcpy(dest->u, src->c, dest->nb_jobs * sizeof(int));
    return val;
}

void partlist_permquicksort(int* perm, PartList* part, int nb_part,
                            int (*functionPtr)(PartList*, PartList*)) {
    int      i, j, temp;
    PartList t;

    if (nb_part <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(nb_part - 1) / 2], temp);
    i = 0;
    j = nb_part;
    memcpy(&t, &(part[perm[0]]), sizeof(PartList));

    while (1) {
        do {
            i++;
        } while (i < nb_part && (*functionPtr)(&(part[perm[i]]), &t));

        do {
            j--;
        } while ((*functionPtr)(&t, &(part[perm[j]])));

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    partlist_permquicksort(perm, part, j, (*functionPtr));
    partlist_permquicksort(perm + i, part, nb_part - i, (*functionPtr));
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
        int        iter = 0;
        // int        Hmax = ((interval*)g_ptr_array_index(intervals,
        // intervals->len - 1))->b;

        for (unsigned i = 0; i < machine->len; ++i) {
            Job*      tmp = (Job*)g_ptr_array_index(machine, i);
            interval* I = (interval*)g_ptr_array_index(intervals, iter);
            while (!(sol->c[tmp->job] <= I->b && I->a < sol->c[tmp->job])) {
                iter++;
                I = (interval*)g_ptr_array_index(intervals, iter);
            }
            sol->u[tmp->job] = iter;
        }
    }
}

void solution_calculate_partition_all(Solution* sol, GPtrArray* intervals) {
    for (int i = 0; i < sol->nb_machines; ++i) {
        solution_calculate_partition_machine(sol, intervals, i);
    }
}

static void calculate_partition(Solution* sol, GPtrArray* intervals, int m,
                                int* u, int* last) {
    int   count = 0;
    int   cur = *last;
    void* tmp;

    GPtrArray* machine = sol->part[m].machine;
    interval*  I = (interval*)g_ptr_array_index(intervals, *u);
    Job*       first_job = g_ptr_array_index(machine, 0);
    Job*       i = (Job*)g_ptr_array_index(machine, cur);
    Job*       j;

    while (sol->c[i->job] > I->a && sol->c[i->job] <= I->b &&
           sol->c[i->job] - i->processing_time <= I->b &&
           sol->c[i->job] - i->processing_time > I->a) {
        if (first_job != i) {
            count++;
            i = (Job*)g_ptr_array_index(machine, --cur);
        } else {
            count++;
            break;
        }
    }

    if (count > 0) {
        if (first_job->job != i->job) {
            g_qsort_with_data(machine->pdata + cur + 1, count, sizeof(Job*),
                              compare_interval, I);
            j = (Job*)g_ptr_array_index(machine, cur + 1);
        } else {
            g_qsort_with_data(machine->pdata, count + 1, sizeof(Job*),
                              compare_interval, I);
            i = (Job*)g_ptr_array_index(machine, 0);
            j = (Job*)g_ptr_array_index(machine, 1);
        }

        if (sol->c[i->job] > I->a && sol->c[i->job] <= I->b) {
            if (compare_interval(&i, &j, I) < 0) {
                cur--;
                *last = cur;
                if (cur >= 0) {
                    i = (Job*)g_ptr_array_index(machine, cur);
                    *u = CC_MIN(*u - 1, sol->u[i->job]);
                } else {
                    *u = -1;
                }
            } else {
                sol->c[j->job] =
                    sol->c[i->job] - i->processing_time + j->processing_time;
                sol->c[i->job] = sol->c[j->job] + i->processing_time;
                CC_SWAP(g_ptr_array_index(machine, cur),
                        g_ptr_array_index(machine, cur + 1), tmp);
                if (sol->c[j->job] > I->a && sol->c[j->job] <= I->b) {
                    g_qsort_with_data(machine->pdata + cur, count + 1,
                                      sizeof(Job*), compare_interval, I);
                    *last = cur;
                    if (cur >= 0) {
                        i = (Job*)g_ptr_array_index(machine, cur);
                        *u = CC_MIN(*u - 1, sol->u[i->job]);
                    } else {
                        *u = -1;
                    }
                }
            }
        }
    } else {
        cur--;
        *last = cur;
        if (cur >= 0) {
            i = (Job*)g_ptr_array_index(machine, cur);
            *u = CC_MIN(*u - 1, sol->u[i->job]);
        } else {
            *u = -1;
        }
    }

    solution_calculate_machine(sol, m);
}

int solution_canonical_order(Solution* sol, GPtrArray* intervals) {
    int val = 0;

    solution_calculate_partition_all(sol, intervals);

    for (int it = 0; it < sol->nb_machines; ++it) {
        GPtrArray* machine = sol->part[it].machine;
        int        last = machine->len - 1;
        Job*       i = (Job*)g_ptr_array_index(machine, last);
        int        u = sol->u[i->job];
        while (u >= 0) {
            calculate_partition(sol, intervals, it, &u, &last);
        }
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
