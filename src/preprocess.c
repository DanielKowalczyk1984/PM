#include <interval.h>
#include <wct.h>

void g_problem_summary_init(gpointer data, gpointer user_data) {
    Job*     j = (Job*)data;
    Problem* prob = (Problem*)user_data;

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

void calculate_Hmax(Problem* problem) {
    int       temp = 0;
    double    temp_dbl = 0.0;
    NodeData* pd = problem->root_pd;

    temp = problem->p_sum - problem->pmax;
    temp_dbl = (double)temp;
    temp_dbl = floor(temp_dbl / problem->nb_machines);
    problem->H_max = pd->H_max = (int)temp_dbl + problem->pmax;
    problem->H_min = pd->H_min =
        (int)ceil(temp_dbl / problem->nb_machines) - problem->pmax;

    GPtrArray* duration = g_ptr_array_copy(problem->g_job_array, NULL, NULL);
    g_ptr_array_set_free_func(duration, NULL);
    g_ptr_array_sort(duration, g_compare_duration);

    int    m = 0;
    int    i = problem->nb_jobs - 1;
    double tmp = problem->p_sum;
    problem->H_min = problem->p_sum;
    do {
        Job* job = g_ptr_array_index(duration, i);
        tmp -= job->processing_time;
        m++;
        i--;

    } while (m < problem->nb_machines - 1);

    problem->H_min = pd->H_min = (int)ceil(tmp / problem->nb_machines);
    g_ptr_array_free(duration, TRUE);
    printf(
        "H_max = %d, H_min = %d,  pmax = %d, pmin = %d, p_sum = %d, off = %d\n",
        problem->H_max, problem->H_min, problem->pmax, problem->pmin,
        problem->p_sum, problem->off);
}

void determine_jobs_order_interval(Problem* problem) {
    interval* tmp_interval;

    GPtrArray* local_intervals = problem->intervals;

    for (unsigned i = 0; i < problem->g_job_array->len; ++i) {
        Job* tmp_j = (Job*)g_ptr_array_index(problem->g_job_array, i);
        // tmp_j->pos_interval = CC_SAFE_MALLOC(local_intervals->len, int);
        for (unsigned j = 0; j < local_intervals->len; ++j) {
            tmp_interval = (interval*)g_ptr_array_index(local_intervals, j);
            GPtrArray* sigma = tmp_interval->sigma;
            for (unsigned k = 0; k < sigma->len; ++k) {
                Job* tmp = (Job*)g_ptr_array_index(sigma, k);
                if (tmp == tmp_j) {
                    // tmp_j->pos_interval[j] = k;
                    break;
                }
            }
        }
    }
}

int preprocess_data(Problem* problem) {
    int       val = 0;
    int       i = 0;
    NodeData* root = problem->root_pd;

    /** Calculate the statistics of the instance */
    g_ptr_array_foreach(problem->g_job_array, g_problem_summary_init, problem);

    /** Calculate H_max */
    calculate_Hmax(problem);

    /** order the jobarray of problem following edd rule */
    g_ptr_array_sort_with_data(problem->g_job_array, g_job_compare_edd, NULL);

    g_ptr_array_foreach(problem->g_job_array, g_set_jobarray_job, &i);
    root->jobarray = problem->g_job_array;
    root->off = problem->off;

    /** Find the intervals of the instance at the root node */
    find_division(problem);

    /** Create all node of the ZDD */
    create_ordered_jobs_array(problem->intervals, root->ordered_jobs);

    /** Determine the position of each job in the interval */
    // determine_jobs_order_interval(problem);

    return val;
}

static int calculate_T(interval_pair* pair, int k, GPtrArray* interval_array) {
    interval* I = (interval*)g_ptr_array_index(interval_array, k);
    Job*      i = pair->a;
    Job*      j = pair->b;
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

static int check_interval(interval_pair* pair,
                          int            k,
                          GPtrArray*     interval_array) {
    interval* I = (interval*)g_ptr_array_index(interval_array, k);
    Job*      j = pair->b;
    return (I->a + j->processing_time >= I->b ||
            calculate_T(pair, k, interval_array) <= I->a);
}

static GPtrArray* array_time_slots(interval* I, GList* pairs) {
    GPtrArray*     array = g_ptr_array_new_with_free_func(free);
    interval_pair* tmp;
    interval_pair* min_data;
    GList*         min;
    int*           tmp_int;

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->a;
    g_ptr_array_add(array, tmp_int);

    while (pairs) {
        min = pairs;
        min_data = (interval_pair*)min->data;
        for (GList* i = min->next; i; i = g_list_next(i)) {
            tmp = ((interval_pair*)i->data);
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
            tmp = (interval_pair*)i->data;
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

void create_ordered_jobs_array(GPtrArray* a, GPtrArray* b) {
    interval*          tmp_interval;
    Job*               tmp_j;
    job_interval_pair* tmp_pair;
    for (unsigned i = 0; i < a->len; ++i) {
        tmp_interval = (interval*)g_ptr_array_index(a, i);
        GPtrArray* jobarray = tmp_interval->sigma;
        for (unsigned j = 0; j < jobarray->len; ++j) {
            tmp_j = (Job*)g_ptr_array_index(jobarray, j);
            if (tmp_j->processing_time <= tmp_interval->b) {
                tmp_pair = CC_SAFE_MALLOC(1, job_interval_pair);
                tmp_pair->j = tmp_j;
                tmp_pair->I = tmp_interval;
                g_ptr_array_add(b, tmp_pair);
            }
        }
    }

    printf("There are %u layers\n", b->len);
}

int find_division(Problem* problem) {
    int            val = 0;
    int            counter = 0;
    int            nb_jobs = problem->nb_jobs;
    int            prev;
    NodeData*      root_pd = problem->root_pd;
    GPtrArray*     tmp_array = g_ptr_array_new_with_free_func(g_interval_free);
    GPtrArray*     jobarray = problem->g_job_array;
    Job*           tmp_j;
    Job *          j1, *j2;
    interval*      tmp_interval;
    interval_pair* pair;
    interval_pair  tmp_pair;

    /** Find initial partition */
    prev = 0;
    for (int i = 0; i < nb_jobs && prev < problem->H_max; ++i) {
        tmp_j = (Job*)g_ptr_array_index(jobarray, i);
        int tmp = CC_MIN(problem->H_max, tmp_j->due_time);
        if (prev < tmp) {
            tmp_interval = interval_alloc(prev, tmp, -1, jobarray, nb_jobs);
            g_ptr_array_add(tmp_array, tmp_interval);
            CCcheck_NULL_2(tmp_interval, "Failed to allocate memory");
            prev = tmp_j->due_time;
        }
    }

    if (prev < problem->H_max) {
        tmp_interval =
            interval_alloc(prev, problem->H_max, -1, jobarray, nb_jobs);
        g_ptr_array_add(tmp_array, tmp_interval);
    }

    /** calculate the new intervals */
    for (unsigned i = 0; i < tmp_array->len; ++i) {
        GList* pairs = (GList*)NULL;
        tmp_interval = (interval*)g_ptr_array_index(tmp_array, i);
        for (size_t j = 0; j < tmp_interval->sigma->len - 1; j++) {
            for (size_t k = j + 1; k < tmp_interval->sigma->len; k++) {
                j1 = (Job*)g_ptr_array_index(tmp_interval->sigma, j);
                j2 = (Job*)g_ptr_array_index(tmp_interval->sigma, k);
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
            GPtrArray* slots;
            slots = array_time_slots(tmp_interval, pairs);
            for (unsigned j = 1; j < slots->len; ++j) {
                g_ptr_array_add(problem->intervals,
                                interval_alloc(*((int*)slots->pdata[j - 1]),
                                               *((int*)slots->pdata[j]),
                                               counter, jobarray, nb_jobs));
                counter++;
            }
            g_ptr_array_free(slots, TRUE);
        } else {
            g_ptr_array_add(problem->intervals,
                            interval_alloc(tmp_interval->a, tmp_interval->b,
                                           counter, jobarray, nb_jobs));
            counter++;
        }
    }

CLEAN:
    g_ptr_array_free(tmp_array, TRUE);
    return val;
}
