#include "wctprivate.h"
#include <fmt/core.h>
#include "interval.h"

static int get_problem_name(char* pname, const char* end_file_name) {
    int           rval = 0;
    unsigned long len = 0;
    const char*   fname = strrchr(end_file_name, '/');
    const char*   lastdot = strrchr(end_file_name, '.');

    if (!fname) {
        fname = end_file_name;
    } else {
        fname++;
    }

    if (lastdot) {
        len = lastdot - fname + 1;
    } else {
        len = strlen(fname);
    }

    if (snprintf(pname, len, "%s", fname) < 0) {
        rval = 1;
    }

    fmt::print("Extracted problem name {}\n", pname);
    return rval;
}

_Problem::~_Problem() { /*free the parameters*/
    parms_free(&(parms));

    g_ptr_array_free(g_job_array, TRUE);
    g_ptr_array_free(ColPool, TRUE);
    g_ptr_array_free(intervals, TRUE);
    solution_free(&(opt_sol));
}

static void g_problem_summary_init(gpointer data, gpointer user_data) {
    Problem* prob = static_cast<Problem*>(user_data);
    Job*     j = static_cast<Job*>(data);

    prob->p_sum += j->processing_time;
    prob->pmax = std::max(prob->pmax, j->processing_time);
    prob->pmin = std::max(prob->pmin, j->processing_time);
    prob->dmax = std::max(prob->dmax, j->due_time);
    prob->dmin = std::max(prob->dmin, j->due_time);
}

void _Problem::calculate_Hmax() {
    int       temp = 0;
    double    temp_dbl = 0.0;
    NodeData* pd = root_pd;

    temp = p_sum - pmax;
    temp_dbl = static_cast<double>(temp);
    temp_dbl = floor(temp_dbl / nb_machines);
    H_max = pd->H_max = static_cast<int>(temp_dbl) + pmax;
    H_min = pd->H_min = static_cast<int>(ceil(temp_dbl / nb_machines)) - pmax;

    GPtrArray* duration = g_ptr_array_copy(g_job_array, NULL, NULL);
    g_ptr_array_set_free_func(duration, NULL);
    g_ptr_array_sort(duration, g_compare_duration);

    int    m = 0;
    int    i = nb_jobs - 1;
    double tmp = p_sum;
    H_min = p_sum;
    do {
        Job* job = static_cast<Job*>(g_ptr_array_index(duration, i));
        tmp -= job->processing_time;
        m++;
        i--;

    } while (m < nb_machines - 1);

    H_min = pd->H_min = static_cast<int>(ceil(tmp / nb_machines));
    g_ptr_array_free(duration, TRUE);
    fmt::print(
        "H_max = {}, H_min = {},  pmax = {}, pmin = {}, p_sum = {}, off = "
        "{}\n",
        H_max, H_min, pmax, pmin, p_sum, off);
}

void _Problem::create_ordered_jobs_array(GPtrArray* a, GPtrArray* b) {
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

int _Problem::find_division() {
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

int _Problem::check_interval(interval_pair* pair,
                             int            k,
                             GPtrArray*     interval_array) {
    auto* I = static_cast<interval*>(g_ptr_array_index(interval_array, k));
    auto* j = pair->b;
    return (I->a + j->processing_time >= I->b ||
            calculate_T(pair, k, interval_array) <= I->a);
}

int _Problem::calculate_T(interval_pair* pair,
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

GPtrArray* _Problem::array_time_slots(interval* I, GList* pairs) {
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

int _Problem::preprocess_data() {
    int       val = 0;
    int       i = 0;
    NodeData* root = root_pd;

    /** Calculate the statistics of the instance */
    g_ptr_array_foreach(g_job_array, g_problem_summary_init, this);

    /** Calculate H_max */
    calculate_Hmax();

    /** order the jobarray of problem following edd rule */
    g_ptr_array_sort_with_data(g_job_array, g_job_compare_edd, NULL);

    g_ptr_array_foreach(g_job_array, g_set_jobarray_job, &i);
    root->jobarray = g_job_array;
    root->off = off;

    /** Find the intervals of the instance at the root node */
    find_division();

    /** Create all node of the ZDD */
    create_ordered_jobs_array(intervals, root_pd->ordered_jobs);

    /** Determine the position of each job in the interval */
    // determine_jobs_order_interval(problem);

    return val;
}
