#include <interval.h>
#include "util.h"

void g_print_interval(gpointer data, MAYBE_UNUSED gpointer user_data) {
    interval* a = (interval*)data;
    printf("interval %d: (%d %d]: ", a->key, a->a, a->b);
    g_ptr_array_foreach(a->sigma, g_print_job, NULL);
    printf("\n");
}

gint g_compare_interval_data(gconstpointer a, gconstpointer b, gpointer data) {
    const Job* x = *(Job* const*)a;
    const Job* y = *(Job* const*)b;
    interval*  user_data = (interval*)data;
    int        diff = user_data->b - user_data->a;
    double     w_x = (x->due_time >= user_data->b)
                         ? 0.0
                         : (double)x->weight / x->processing_time;
    double     w_y = (y->due_time >= user_data->b)
                         ? 0.0
                         : (double)y->weight / y->processing_time;

    if (x->processing_time >= diff) {
        if (y->processing_time < diff) {
            return -1;
        } else {
            if (w_x > w_y) {
                return -1;
            } else if (w_y > w_x) {
                return 1;
            } else if (x->processing_time > y->processing_time) {
                return -1;
            } else if (y->processing_time > x->processing_time) {
                return 1;
            }

            return x->job - y->job;
        }
    } else {
        if (y->processing_time >= diff) {
            return 1;
        } else {
            if (w_x > w_y) {
                return -1;
            } else if (w_y > w_x) {
                return 1;
            } else if (x->processing_time > y->processing_time) {
                return -1;
            } else if (y->processing_time > x->processing_time) {
                return 1;
            }

            return x->job - y->job;
        }
    }
}

void interval_init(interval*  p,
                   int        a,
                   int        b,
                   int        key,
                   GPtrArray* jobarray,
                   int        nb_jobs) {
    p->a = a;
    p->b = b;
    p->key = key;
    p->sigma = g_ptr_array_copy(jobarray, NULL, NULL);
    g_ptr_array_set_free_func(p->sigma, NULL);
    p->begin = 0;

    g_ptr_array_sort_with_data(p->sigma, g_compare_interval_data, p);

    Job* j = (Job*)g_ptr_array_index(p->sigma, p->begin);
    while (p->b - p->a <= j->processing_time && p->begin < (int)jobarray->len) {
        j = (Job*)g_ptr_array_index(p->sigma, p->begin);
        p->begin++;
    }
}

interval* interval_alloc(int        a,
                         int        b,
                         int        key,
                         GPtrArray* jobarray,
                         int        nb_jobs) {
    interval* p = CC_SAFE_MALLOC(1, interval);
    CCcheck_NULL_3(p, "Failed to allocate memory");
    interval_init(p, a, b, key, jobarray, nb_jobs);
CLEAN:
    return p;
}

// interval* interval_copy(interval* src) {
//     interval* ret = CC_SAFE_MALLOC(1, interval);
//     CCcheck_NULL_3(ret, "Failed to allocate memory");

//     ret->a = src->a;
//     ret->b = src->b;
//     ret->begin = src->begin;

//     g_ptr_array_sized_new(src->sigma->len);

//     for (unsigned i = 0; i < src->sigma->len; ++i) {
//         g_ptr_array_add(ret->sigma, g_ptr_array_index(src->sigma, i));
//     }
// CLEAN:
//     return ret;
// }

gpointer g_copy_interval(gconstpointer src, gpointer data) {
    interval*       aux = CC_SAFE_MALLOC(1, interval);
    const interval* src_interval = (const interval*)src;

    aux->a = src_interval->a;
    aux->b = src_interval->b;
    aux->begin = src_interval->begin;

    aux->sigma = g_ptr_array_copy(src_interval->sigma, NULL, NULL);

    return aux;
}

gpointer g_copy_interval_pair(gconstpointer src, gpointer data) {
    job_interval_pair*       aux = CC_SAFE_MALLOC(1, job_interval_pair);
    const job_interval_pair* pair = (const job_interval_pair*)src;

    aux->I = pair->I;
    aux->j = pair->j;

    return aux;
}

void interval_free(interval* p) {
    if (p != (interval*)NULL) {
        g_ptr_array_free(p->sigma, TRUE);
    }
}

void g_interval_free(void* p) {
    if (p != NULL) {
        interval* part = (interval*)p;
        g_ptr_array_free(part->sigma, TRUE);
        CC_FREE(part, interval);
    }
}

void interval_pair_free(void* p) {
    if (p != NULL) {
        interval_pair* tmp = (interval_pair*)p;
        CC_FREE(tmp, interval_pair);
    }
}

void print_interval_pair(GPtrArray* ordered_jobs) {
    interval*          cur = (interval*)NULL;
    job_interval_pair* tmp_p = NULL;

    for (size_t i = 0; i < ordered_jobs->len; i++) {
        tmp_p =
            static_cast<job_interval_pair*>(g_ptr_array_index(ordered_jobs, i));
        if (tmp_p->I != cur) {
            cur = tmp_p->I;
            printf("\n");
            printf("Interval %d (%d,%d]: ", cur->key, cur->a, cur->b);
        }
        printf("%d ", tmp_p->j->job);
    }
    printf("\n");
}

void count_jobs_interval_pair(GPtrArray* ordered_jobs) {
    job_interval_pair* tmp_p = NULL;

    for (size_t i = 0; i < ordered_jobs->len; i++) {
        tmp_p =
            static_cast<job_interval_pair*>(g_ptr_array_index(ordered_jobs, i));
        tmp_p->j->num_layers++;
    }
}
