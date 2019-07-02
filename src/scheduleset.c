
////////////////////////////////////////////////////////////////
//                                                            //
//  scheduleset.c                                                //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <defs.h>
#include <scheduleset.h>
#include <util.h>
void iterator(gpointer key, gpointer value, gpointer user_data);

void iterator(gpointer key, gpointer value, gpointer user_data) {
    GHashTable *new_table = (GHashTable *)user_data;
    g_hash_table_insert(new_table, key, value);
}

void scheduleset_init(ScheduleSet *set) {
    if (set) {
        set->num = (int *)NULL;
        set->del = 0;
        set->age = 0;
        set->total_processing_time = 0;
        set->total_weighted_completion_time = 0;
        set->id = -1;
        set->total_weighted_completion_time = 0;
        set->job_list = g_ptr_array_new();
        set->edge_list = g_ptr_array_new();
        set->table = g_hash_table_new(g_direct_hash, g_direct_equal);
    }
}

void scheduleset_init_bis(ScheduleSet *set) {
    if (set) {
        set->num = (int *)NULL;
        set->del = 0;
        set->age = 0;
        set->total_processing_time = 0;
        set->total_weighted_completion_time = 0;
        set->id = -1;
        set->total_weighted_completion_time = 0;
        set->job_list = NULL;
        set->edge_list = NULL;
        set->table = g_hash_table_new(g_direct_hash, g_direct_equal);
    }
}

void scheduleset_free(ScheduleSet *set) {
    if (set) {
        CC_IFFREE(set->num, int);
        if(set->job_list) {
            g_ptr_array_free(set->job_list, TRUE);
        }
        if(set->edge_list) {
            g_ptr_array_free(set->edge_list, TRUE);
        }
        g_hash_table_destroy(set->table);

        set->total_processing_time = 0;
        set->age = 0;
        set->del = 0;
        set->total_weighted_completion_time = 0;
        CC_IFFREE(set->num, int);
    }
}

void g_scheduleset_free(void *set) {
    ScheduleSet *tmp = (ScheduleSet *)set;
    if (tmp) {
        CC_IFFREE(tmp->num, int);

        tmp->total_processing_time = 0;
        tmp->age = 0;
        tmp->del = 0;
        tmp->total_weighted_completion_time = 0;
        tmp->id = -1;
        tmp->total_weighted_completion_time = 0;
        if(tmp->job_list) {
            g_ptr_array_free(tmp->job_list, TRUE);
        }
        if(tmp->edge_list) {
            g_ptr_array_free(tmp->edge_list, TRUE);
        }
        g_hash_table_destroy(tmp->table);
        CC_IFFREE(tmp, ScheduleSet);
    }
}

ScheduleSet *scheduleset_alloc(int nbjobs) {
    ScheduleSet *tmp;
    tmp = CC_SAFE_MALLOC(1, ScheduleSet);
    CCcheck_NULL_3(tmp, "Failed to allocate memory") scheduleset_init(tmp);
    tmp->num = CC_SAFE_MALLOC(nbjobs, int);
    fill_int(tmp->num, nbjobs, 0);

CLEAN:
    return tmp;
}

ScheduleSet *scheduleset_alloc_bis(int nbjobs) {
    ScheduleSet *tmp;
    tmp = CC_SAFE_MALLOC(1, ScheduleSet);
    CCcheck_NULL_3(tmp, "Failed to allocate memory") scheduleset_init_bis(tmp);
    tmp->num = CC_SAFE_MALLOC(nbjobs, int);
    fill_int(tmp->num, nbjobs, 0);

CLEAN:
    return tmp;
}

void g_sum_processing_time(gpointer data, gpointer user_data) {
    Job *        j = (Job *)data;
    ScheduleSet *set = (ScheduleSet *)user_data;

    set->total_processing_time += j->processing_time;
    set->total_weighted_completion_time += value_Fj(set->total_processing_time, j);
    (set->num[j->job])++;
    g_ptr_array_add(set->job_list, j);
}

ScheduleSet *scheduleset_from_solution(GPtrArray *machine, int nbjobs) {
    ScheduleSet *tmp;

    tmp = CC_SAFE_MALLOC(1, ScheduleSet);
    CCcheck_NULL_3(tmp, "failed to allocate memory")

    scheduleset_init(tmp);
    tmp->num = CC_SAFE_MALLOC(nbjobs, int);
    CCcheck_NULL(tmp->num, "Failed to allocate memory")
    fill_int(tmp->num, nbjobs, 0);
    g_ptr_array_foreach(machine, g_sum_processing_time, tmp);

CLEAN:
    return tmp;
}

void schedulesets_free(ScheduleSet **sets, int *nsets) {
    if (*sets) {
        for (int i = 0; i < *nsets; i++) {
            scheduleset_free(&(*sets)[i]);
        }

        CC_IFFREE(*sets, ScheduleSet);
    }

    *nsets = 0;
}

void scheduleset_SWAP(ScheduleSet *c1, ScheduleSet *c2, ScheduleSet *t) {
    if (c1 != c2) {
        memcpy(t, c2, sizeof(ScheduleSet));
        memcpy(c2, c1, sizeof(ScheduleSet));
        memcpy(c1, t, sizeof(ScheduleSet));
    }
}

int scheduleset_less(ScheduleSet *c1, ScheduleSet *c2) {
    int        i;
    GPtrArray *tmp1 = c1->job_list;
    GPtrArray *tmp2 = c2->job_list;

    if (tmp1->len != tmp2->len) {
        return tmp1->len - tmp2->len;
    }

    for (i = 0; i < tmp1->len; ++i) {
        Job* tmp_j1 = (Job *)g_ptr_array_index(tmp1, i);
        Job* tmp_j2 = (Job *)g_ptr_array_index(tmp2, i);
        if (tmp_j1->job != tmp_j2->job) {
            return tmp_j1->job - tmp_j2->job;
        }
    }

    return 0;
}

gint g_scheduleset_less(gconstpointer a, gconstpointer b) {
    int                i;
    const ScheduleSet *c1 = *((ScheduleSet *const *)a);
    const ScheduleSet *c2 = *((ScheduleSet *const *)b);
    GPtrArray *        tmp1 = c1->job_list;
    GPtrArray *        tmp2 = c2->job_list;

    if (tmp1->len != tmp2->len) {
        return tmp1->len - tmp2->len;
    }

    for (i = 0; i < tmp1->len; ++i) {
        Job* tmp_j1 = (Job *)g_ptr_array_index(tmp1, i);
        Job* tmp_j2 = (Job *)g_ptr_array_index(tmp2, i);
        if (tmp_j1->job != tmp_j2->job) {
            return tmp_j1->job - tmp_j2->job;
        }
    }

    return 0;
}

void g_scheduleset_print(gpointer data, gpointer user_data) {
    ScheduleSet *tmp = (ScheduleSet *)data;
    GPtrArray *  tmp_a = tmp->job_list;
    printf("Machine %d: ", tmp->id);

    g_ptr_array_foreach(tmp_a, g_print_machine, NULL);

    printf("with C = %d, cost = %d and %u jobs\n", tmp->total_processing_time, tmp->total_weighted_completion_time,
           tmp_a->len);
}

void g_compute_nblayers_schedule(gpointer data, gpointer user_data){
    Job *j = (Job *) data;
    ScheduleSet *tmp = (ScheduleSet *) user_data;
    if(tmp->num[j->job] > 1) {
        j->num_layers = 1;
    }
}

int print_schedule(ScheduleSet *cclasses, int ccount) {
    int i;
    int sum = 0;

    for (i = 0; i < ccount; i++) {
        printf("Machine %d:", i);

        g_ptr_array_foreach(cclasses[i].job_list, g_print_machine, NULL);

        printf(" with total_processing_time %d and %d jobs\n", cclasses[i].total_processing_time,
               cclasses[i].job_list->len);
        sum += cclasses[i].job_list->len;
    }

    printf("Total of jobs = %d\n", sum);
    fflush(stdout);
    fflush(stdout);
    return 0;
}

int scheduleset_max(ScheduleSet *cclasses, int ccount) {
    int val = 0;
    int i;

    for (i = 0; i < ccount; i++) {
        if (cclasses[i].total_processing_time > val) {
            val = cclasses[i].total_processing_time;
        }
    }

    return val;
}
