#include "job.h"
#include "wctprivate.h"

void g_problem_summary_init(gpointer data, gpointer user_data) {
    Job*     j = static_cast<Job*>(data);
    Problem* prob = static_cast<Problem*>(user_data);

    prob->p_sum += j->processing_time;
    prob->pmax = CC_MAX(prob->pmax, j->processing_time);
    prob->pmin = CC_MIN(prob->pmin, j->processing_time);
    prob->dmax = CC_MAX(prob->dmax, j->due_time);
    prob->dmin = CC_MIN(prob->dmin, j->due_time);
}

gint g_job_compare_edd(const void* a, const void* b, MAYBE_UNUSED void* data) {
    const Job* x = *((Job* const*)a);
    const Job* y = *((Job* const*)b);

    if (x->due_time > y->due_time) {
        return (1);
    } else if (x->due_time < y->due_time) {
        return (-1);
    } else if (x->processing_time > y->processing_time) {
        return (1);
    } else if (x->processing_time < y->processing_time) {
        return (-1);
    } else if (x->weight < y->weight) {
        return (1);
    } else if (x->weight > y->weight) {
        return (-1);
    } else if (x->job > y->job) {
        return (1);
    } else if (x->job < y->job) {
        return (-1);
    }

    return (0);
}

gint g_compare_duration(gconstpointer a, gconstpointer b) {
    const Job* x = *((Job* const*)a);
    const Job* y = *((Job* const*)b);

    if (x->processing_time < y->processing_time) {
        return -1;
    } else {
        return 1;
    }
}

int Problem::calculate_T(interval_pair* pair,
                         int            k,
                         GPtrArray*     interval_array) {
    auto* I = static_cast<interval*>(g_ptr_array_index(interval_array, k));
    auto* i = pair->a;
    auto* j = pair->b;
    pair->left = I->a;
    pair->right = I->a + j->processing_time;

    if (pair->left > I->b - i->processing_time) {
        return pair->left;
    } else {
        if (value_diff_Fij(pair->left, i, j) <= 0) {
            return pair->left;
        }

        // for (int t = k - 1; t >= 0; t--) {
        //     tmp = (interval*)g_ptr_array_index(interval_array, t);
        //     pair->left = tmp->a + j->processing_time - i->processing_time;

        //     if (value_diff_Fij(pair->left, i, j) <= 0 && pair->left >= tmp->a
        //     &&
        //         pair->left <= tmp->b - i->processing_time) {
        //         break;
        //     }
        // }

        pair->left =
            i->due_time +
            (int)ceil((double)(j->weight * i->processing_time) / i->weight) -
            i->processing_time;
        return pair->left;
    }
}

int Problem::check_interval(interval_pair* pair,
                            int            k,
                            GPtrArray*     interval_array) {
    auto* I = static_cast<interval*>(g_ptr_array_index(interval_array, k));
    auto* j = pair->b;
    return (I->a + j->processing_time >= I->b ||
            calculate_T(pair, k, interval_array) <= I->a);
}

GPtrArray* Problem::array_time_slots(interval* I, GList* pairs) {
    GPtrArray* array = g_ptr_array_new_with_free_func(free);
    // interval_pair* tmp = (interval_pair*)NULL;
    interval_pair* min_data = nullptr;
    GList*         min = (GList*)NULL;
    int*           tmp_int = NULL;

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->a;
    g_ptr_array_add(array, tmp_int);

    while (pairs) {
        min = pairs;
        min_data = static_cast<interval_pair*>(min->data);
        for (GList* i = min->next; i; i = g_list_next(i)) {
            auto* tmp = static_cast<interval_pair*>(i->data);
            if (tmp->right < min_data->right) {
                min = i;
                min_data = (interval_pair*)i->data;
            }
        }

        tmp_int = CC_SAFE_MALLOC(1, int);
        *tmp_int = min_data->right;
        g_ptr_array_add(array, tmp_int);
        pairs = g_list_remove_link(pairs, min);
        g_list_free_full(min, interval_pair_free);

        GList* i = pairs;
        while (i) {
            auto* tmp = static_cast<interval_pair*>(i->data);
            if ((*tmp_int >= tmp->left && *tmp_int <= tmp->right) ||
                (*tmp_int + tmp->b->processing_time >= I->b)) {
                GList* remove = i;
                i = g_list_next(i);
                pairs = g_list_remove_link(pairs, remove);
                g_list_free_full(remove, interval_pair_free);
            } else {
                tmp->right = *tmp_int + tmp->b->processing_time;
                i = g_list_next(i);
            }
        }
    }

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->b;
    g_ptr_array_add(array, tmp_int);

    return array;
}

int Problem::find_division() {
    int            val = 0;
    int            counter = 0;
    int            prev = 0;
    GPtrArray*     tmp_array = g_ptr_array_new_with_free_func(g_interval_free);
    Job*           tmp_j = NULL;
    Job *          j1 = NULL, *j2 = NULL;
    interval*      tmp_interval = NULL;
    interval_pair* pair = NULL;
    interval_pair  tmp_pair;

    /** Find initial partition */
    prev = 0;
    for (int i = 0; i < nb_jobs && prev < H_max; ++i) {
        tmp_j = static_cast<Job*>(g_ptr_array_index(g_job_array, i));
        int tmp = CC_MIN(H_max, tmp_j->due_time);
        if (prev < tmp) {
            tmp_interval = interval_alloc(prev, tmp, -1, g_job_array, nb_jobs);
            g_ptr_array_add(tmp_array, tmp_interval);
            CCcheck_NULL_2(tmp_interval, "Failed to allocate memory");
            prev = tmp_j->due_time;
        }
    }

    if (prev < H_max) {
        tmp_interval = interval_alloc(prev, H_max, -1, g_job_array, nb_jobs);
        g_ptr_array_add(tmp_array, tmp_interval);
    }

    /** calculate the new intervals */
    for (unsigned i = 0; i < tmp_array->len; ++i) {
        GList* pairs = (GList*)NULL;
        tmp_interval = static_cast<interval*>(g_ptr_array_index(tmp_array, i));
        for (size_t j = 0; j < tmp_interval->sigma->len - 1; j++) {
            for (size_t k = j + 1; k < tmp_interval->sigma->len; k++) {
                j1 = static_cast<Job*>(
                    g_ptr_array_index(tmp_interval->sigma, j));
                j2 = static_cast<Job*>(
                    g_ptr_array_index(tmp_interval->sigma, k));
                tmp_pair =
                    (interval_pair){j1, j2, tmp_interval->a, tmp_interval->b};
                if (!check_interval(&tmp_pair, i, tmp_array) &&
                    !(j1->due_time >= tmp_interval->b) &&
                    !(j2->due_time >= tmp_interval->b)) {
                    pair = CC_SAFE_MALLOC(1, interval_pair);
                    *pair =
                        (interval_pair){j1, j2, tmp_pair.left, tmp_pair.right};
                    pairs = g_list_append(pairs, pair);
                }
            }
        }

        if (pairs) {
            GPtrArray* slots = array_time_slots(tmp_interval, pairs);
            for (unsigned j = 1; j < slots->len; ++j) {
                g_ptr_array_add(
                    intervals, interval_alloc(*((int*)slots->pdata[j - 1]),
                                              *((int*)slots->pdata[j]), counter,
                                              g_job_array, nb_jobs));
                counter++;
            }
            g_ptr_array_free(slots, TRUE);
        } else {
            g_ptr_array_add(intervals,
                            interval_alloc(tmp_interval->a, tmp_interval->b,
                                           counter, g_job_array, nb_jobs));
            counter++;
        }
    }

CLEAN:
    g_ptr_array_free(tmp_array, TRUE);
    return val;
}

void Problem::create_ordered_jobs_array(GPtrArray* a, GPtrArray* b) {
    // interval* tmp_interval = (interval*)NULL;
    // Job*               tmp_j = nullptr;
    // job_interval_pair* tmp_pair = static_cast<job_interval_pair*)NULL;
    for (unsigned i = 0; i < a->len; ++i) {
        auto* tmp_interval = static_cast<interval*>(g_ptr_array_index(a, i));
        GPtrArray* jobarray = tmp_interval->sigma;
        for (unsigned j = 0; j < jobarray->len; ++j) {
            auto* tmp_j = static_cast<Job*>(g_ptr_array_index(jobarray, j));
            if (tmp_j->processing_time <= tmp_interval->b) {
                job_interval_pair* tmp_pair =
                    CC_SAFE_MALLOC(1, job_interval_pair);
                tmp_pair->j = tmp_j;
                tmp_pair->I = tmp_interval;
                g_ptr_array_add(b, tmp_pair);
            }
        }
    }

    fmt::print("There are {} layers\n", b->len);
}

int Problem::preprocess_data() {
    int val = 0;
    int i = 0;

    /** Calculate the statistics of the instance */
    g_ptr_array_foreach(g_job_array, g_problem_summary_init, this);

    /** Calculate H_max */
    calculate_Hmax();

    /** order the jobarray of problem following edd rule */
    g_ptr_array_sort_with_data(g_job_array, g_job_compare_edd, NULL);

    g_ptr_array_foreach(g_job_array, g_set_jobarray_job, &i);
    root_pd->jobarray = g_job_array;
    root_pd->off = off;

    /** Find the intervals of the instance at the root node */
    find_division();

    /** Create all node of the ZDD */
    create_ordered_jobs_array(intervals, root_pd->ordered_jobs);

    /** Determine the position of each job in the interval */
    // determine_jobs_order_interval(problem);

    return val;
}
