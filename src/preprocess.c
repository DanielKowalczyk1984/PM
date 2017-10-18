#include <interval.h>
#include <wct.h>

void g_problem_summary_init(gpointer data, gpointer user_data) {
    Job *       j = (Job *)data;
    wctproblem *prob = (wctproblem *)user_data;

    prob->psum += j->processingime;
    prob->pmax = CC_MAX(prob->pmax, j->processingime);
    prob->pmin = CC_MIN(prob->pmin, j->processingime);
    prob->dmax = CC_MAX(prob->dmax, j->duetime);
    prob->dmin = CC_MIN(prob->pmin, j->duetime);
}

gint g_job_compare_edd(const void *a, const void *b, void *data) {
    const Job *x = *((Job *const *)a);
    const Job *y = *((Job *const *)b);

    if (x->duetime > y->duetime) {
        return (1);
    } else if (x->duetime < y->duetime) {
        return (-1);
    } else if (x->processingime > y->processingime) {
        return (1);
    } else if (x->processingime < y->processingime) {
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

int calculate_Hmin(
    int *durations, int nmachines, int njobs, int *perm, double *H) {
    int    i, val = 0;
    double temp;

    for (i = 0; i < njobs; ++i) {
        val += durations[i];
    }

    for (i = 0; i < nmachines - 1; ++i) {
        val -= durations[perm[i]];
    }

    temp = (double)val;
    *H = temp / (double)nmachines;
    val = (int)floor(*H);
    return val;
}

void calculate_Hmax(wctproblem *problem) {
    int      temp = 0;
    double   temp_dbl = 0.0;
    wctdata *pd = &(problem->root_pd);

    temp = problem->psum - problem->pmax;
    temp_dbl = (double)temp;
    temp_dbl = floor(temp_dbl / problem->nmachines);
    problem->H_max = pd->H_max = (int)temp_dbl + problem->pmax;
    problem->H_min = (int)ceil(temp_dbl / problem->nmachines) - problem->pmax;
    printf("H_max = %d,  pmax = %d, pmin = %d, psum = %d, off = %d\n",
           problem->H_max, problem->pmax, problem->pmin, problem->psum,
           problem->off);
}

void determine_jobs_order_interval(wctproblem *problem) {
    Job *     tmp, *tmp_j;
    interval *tmp_interval;

    GPtrArray *local_intervals = problem->root_pd.local_intervals;

    for (unsigned i = 0; i < problem->g_job_array->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(problem->g_job_array, i);
        tmp_j->pos_interval = CC_SAFE_MALLOC(local_intervals->len, int);
        for (unsigned j = 0; j < local_intervals->len; ++j) {
            tmp_interval = (interval *)g_ptr_array_index(local_intervals, j);
            GPtrArray *sigma = tmp_interval->sigma;
            for (unsigned k = 0; k < sigma->len; ++k) {
                tmp = (Job *)g_ptr_array_index(sigma, k);
                if (tmp == tmp_j) {
                    tmp_j->pos_interval[j] = k;
                    break;
                }
            }
        }
    }
}

int preprocess_data(wctproblem *problem) {
    int val = 0;
    int i = 0;

    /** Calculate the statistics of the instance */
    g_ptr_array_foreach(problem->g_job_array, g_problem_summary_init, problem);

    /** Calculate H_max */
    calculate_Hmax(problem);

    /** order the jobarray of problem following edd rule */
    g_ptr_array_sort_with_data(problem->g_job_array, g_job_compare_edd, NULL);

    g_ptr_array_foreach(problem->g_job_array, g_set_jobarray_job, &i);

    /** Find the intervals of the instance at the root node */
    find_division(problem);

    /** Create all node of the ZDD */
    create_ordered_jobs_array(problem->root_pd.local_intervals,
                              problem->root_pd.ordered_jobs);

    /** Determine the position of each job in the interval */
    determine_jobs_order_interval(problem);

    return val;
}

static int calculate_T(interval_pair *pair, int k, GPtrArray *interval_array) {
    interval *I = (interval *)g_ptr_array_index(interval_array, k);
    interval *tmp;
    Job *     i = pair->a;
    Job *     j = pair->b;
    pair->left = I->a;
    pair->right = I->a + j->processingime;

    if (pair->left > I->b - i->processingime) {
        return pair->left;
    } else {
        if (value_diff_Fij(pair->left, i, j) <= 0) {
            return pair->left;
        }

        for (int t = k - 1; t >= 0; t--) {
            tmp = (interval *)g_ptr_array_index(interval_array, t);
            pair->left = tmp->a + j->processingime - i->processingime;

            if (value_diff_Fij(pair->left, i, j) <= 0) {
                break;
            }
        }

        return pair->left;
    }
}

static int check_interval(interval_pair *pair,
                          int            k,
                          GPtrArray *    interval_array) {
    interval *I = (interval *)g_ptr_array_index(interval_array, k);
    Job *     j = pair->b;
    return (I->a + j->processingime >= I->b ||
            calculate_T(pair, k, interval_array) <= I->a);
}

static GPtrArray *array_time_slots(interval *I, GList *pairs) {
    GPtrArray *    array = g_ptr_array_new_with_free_func(free);
    interval_pair *tmp;
    interval_pair *min_data;
    GList *        min;
    int *          tmp_int, prev;

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->a;
    g_ptr_array_add(array, tmp_int);
    prev = *tmp_int;

    while (pairs) {
        min = pairs;
        min_data = (interval_pair *)min->data;
        for (GList *i = min->next; i; i = g_list_next(i)) {
            tmp = ((interval_pair *)i->data);
            if (tmp->right < min_data->right) {
                min = i;
                min_data = (interval_pair *)i->data;
            }
        }

        tmp_int = CC_SAFE_MALLOC(1, int);
        *tmp_int = min_data->right;
        g_ptr_array_add(array, tmp_int);
        pairs = g_list_remove_link(pairs, min);
        g_list_free_full(min, interval_pair_free);

        GList *i = pairs;
        while (i) {
            tmp = (interval_pair *)i->data;
            if (*tmp_int >= tmp->left && *tmp_int <= tmp->right) {
                GList *remove = i;
                i = g_list_next(i);
                pairs = g_list_remove_link(pairs, remove);
                g_list_free_full(remove, interval_pair_free);
            } else {
                tmp->right += *tmp_int - prev;
                i = g_list_next(i);
            }
        }

        prev = *tmp_int;
    }

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->b;
    g_ptr_array_add(array, tmp_int);

    return array;
}

void create_ordered_jobs_array(GPtrArray *a, GPtrArray *b) {
    interval *         tmp_interval;
    Job *              tmp_j;
    GPtrArray *        jobarray;
    job_interval_pair *tmp_pair;
    for (unsigned i = 0; i < a->len; ++i) {
        tmp_interval = (interval *)g_ptr_array_index(a, i);
        jobarray = tmp_interval->sigma;
        for (unsigned j = 0; j < jobarray->len; ++j) {
            tmp_j = (Job *)g_ptr_array_index(jobarray, j);
            if (tmp_j->processingime <= tmp_interval->b) {
                tmp_pair = CC_SAFE_MALLOC(1, job_interval_pair);
                tmp_pair->j = tmp_j;
                tmp_pair->I = tmp_interval;
                g_ptr_array_add(b, tmp_pair);
            }
        }
    }
}

int find_division(wctproblem *problem) {
    int            val = 0;
    int            counter = 0;
    int            njobs = problem->njobs;
    int            tmp;
    int            prev;
    wctdata *      root_pd = &(problem->root_pd);
    GPtrArray *    tmp_array = g_ptr_array_new_with_free_func(g_interval_free);
    GPtrArray *    jobarray = problem->g_job_array;
    Job *          tmp_j;
    Job *          j1, *j2;
    interval *     tmp_interval;
    interval_pair *pair;
    interval_pair  tmp_pair;

    /** Find initial partition */
    prev = 0;
    for (int i = 0; i < njobs && prev < problem->H_max; ++i) {
        tmp_j = (Job *)g_ptr_array_index(jobarray, i);
        tmp = CC_MIN(problem->H_max, tmp_j->duetime);
        if (prev < tmp) {
            tmp_interval = interval_alloc(prev, tmp, -1, jobarray, njobs);
            g_ptr_array_add(tmp_array, tmp_interval);
            CCcheck_NULL_2(tmp_interval, "Failed to allocate memory");
            prev = tmp_j->duetime;
        }
    }

    if (prev < problem->H_max) {
        tmp_interval =
            interval_alloc(prev, problem->H_max, -1, jobarray, njobs);
        g_ptr_array_add(tmp_array, tmp_interval);
    }

    /** calculate the new intervals */
    for (unsigned i = 0; i < tmp_array->len; ++i) {
        GList *pairs = (GList *)NULL;
        tmp_interval = (interval *)g_ptr_array_index(tmp_array, i);
        for (size_t j = 0; j < tmp_interval->sigma->len - 1; j++) {
            for (size_t k = j + 1; k < tmp_interval->sigma->len; k++) {
                j1 = (Job *)g_ptr_array_index(tmp_interval->sigma, j);
                j2 = (Job *)g_ptr_array_index(tmp_interval->sigma, k);
                tmp_pair = (interval_pair){j1, j2};
                if (!check_interval(&tmp_pair, i, tmp_array)) {
                    pair = CC_SAFE_MALLOC(1, interval_pair);
                    *pair =
                        (interval_pair){j1, j2, tmp_pair.left, tmp_pair.right};
                    pairs = g_list_append(pairs, pair);
                }
            }
        }

        if (pairs) {
            GPtrArray *slots;
            slots = array_time_slots(tmp_interval, pairs);
            for (unsigned j = 1; j < slots->len; ++j) {
                g_ptr_array_add(root_pd->local_intervals,
                                interval_alloc(*((int *)slots->pdata[j - 1]),
                                               *((int *)slots->pdata[j]),
                                               counter, jobarray, njobs));
                counter++;
            }
            g_ptr_array_free(slots, TRUE);
        } else {
            g_ptr_array_add(root_pd->local_intervals,
                            interval_alloc(tmp_interval->a, tmp_interval->b,
                                           counter, jobarray, njobs));
            counter++;
        }
    }

CLEAN:
    g_ptr_array_free(tmp_array, TRUE);
    return val;
}
