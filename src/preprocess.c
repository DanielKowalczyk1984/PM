#include <wct.h>
#include <interval.h>

void g_problem_summary_init(gpointer data, gpointer user_data){
    Job *j = (Job *) data;
    wctproblem *prob = (wctproblem*) user_data;

    prob->psum += j->processingime;
    prob->pmax = CC_MAX(prob->pmax, j->processingime);
    prob->pmin = CC_MIN(prob->pmin, j->processingime);
    prob->dmax = CC_MAX(prob->dmax, j->duetime);
    prob->dmin = CC_MIN(prob->pmin, j->duetime);
}

gint g_job_compare_edd(const void *a, const void *b, void *data) {
    const Job *x = *((Job * const *)a);
    const Job *y = *((Job * const *)b);

    if (x->duetime > y->duetime) {
        return (1);
    } else if (x->duetime < y->duetime) {
        return (-1);
    } else if (x->processingime > y->processingime) {
        return (1);
    } else if (x->processingime < y->processingime) {
        return (-1);
    } else if (x->weight > y->weight) {
        return (1);
    } else if (x->weight < y->weight) {
        return (-1);
    } else if (x->job > y->job) {
        return (1);
    } else if (x->job < y->job) {
        return (-1);
    }

    return (0);
}

int calculate_Hmax(Job *jobarray, int nmachines, int njobs) {
    int    i, max = jobarray[0].processingime, val = 0;
    double temp;

    for (i = 0; i < njobs; ++i) {
        max = CC_MAX(jobarray[i].processingime, max);
        val += jobarray[i].processingime;
    }

    val -= max;
    temp = (double)val;
    temp = temp / (double)nmachines + max;
    val = (int)ceil(temp);
    return val;
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

int preprocess_data(wctproblem *problem) {
    int      val = 0;
    int      temp = 0;
    double   temp_dbl = 0.0;
    wctdata *pd = &(problem->root_pd);

    g_ptr_array_foreach(problem->g_job_array, g_problem_summary_init, problem);

    /** Calculate H_max */
    temp = problem->psum - problem->pmax;
    temp_dbl = (double)temp;
    temp_dbl = temp_dbl / problem->nmachines + problem->pmax;
    problem->H_max = pd->H_max = (int)ceil(temp_dbl);
    printf("H_max = %d\n", problem->H_max);

    g_ptr_array_sort_with_data(problem->g_job_array, g_job_compare_edd, NULL);

    find_division(problem);

    return val;
}

static int calculate_T(Job *i, Job *j, int k, GPtrArray *interval_array){
    interval *I = (interval *) g_ptr_array_index(interval_array, k);

    int T = I->a;
    if (T > I->b - i->processingime) {
        return T;
    } else {
        T = I->a;

        if(value_diff_Fij(T, i, j) <= 0) {
            return T;
        }
        for(int t = k - 1; t >= 0 ; t--) {
            T = ((interval *) g_ptr_array_index(interval_array, t))->a + j->processingime - i->processingime;
            if(value_diff_Fij(T, i, j) <= 0 || T == 0) {
                break;
            }
        }

        return T;
    }
}

static int check_interval(Job *i, Job *j, int k,GPtrArray *interval_array){
    interval *I = (interval *) g_ptr_array_index(interval_array, k);
    return (I->a + j->processingime >= I->b || calculate_T(i, j, k, interval_array) <= I->a);
}

static int check_interval2(Job *i, Job *j, int k, GPtrArray *interval_array){
    interval *I = (interval *) g_ptr_array_index(interval_array, k);
    printf("test %d %d %d %d\n",I->a, calculate_T(i, j, k, interval_array),I->a + j->processingime, I->b);
    printf("jobs %f %f\n",(double) i->weight/i->processingime, (double) j->weight/j->processingime );
    return (I->a + j->processingime <= I->b && calculate_T(i, j, k, interval_array) >= I->a);
}

int find_division(wctproblem *problem){
    int val = 0;
    int njobs = problem->njobs;
    int tmp;
    int prev;
    GPtrArray *tmp_array = g_ptr_array_new_with_free_func(intervals_free);
    GPtrArray* jobarray = problem->g_job_array;
    Job *tmp_j;
    Job *j1,*j2;
    interval *tmp_interval;

    /** Find initial partition */
    prev = 0;
    for(int i = 0; i < njobs && prev < problem->H_max; ++i) {
        tmp_j = (Job *)g_ptr_array_index(jobarray, i);
        tmp = CC_MIN(problem->H_max, tmp_j->duetime);
        if(prev < tmp) {
            tmp_interval = interval_alloc(prev, tmp, jobarray, njobs);
            g_ptr_array_add(tmp_array, tmp_interval);
            CCcheck_NULL_2(tmp_interval, "Failed to allocate memory")
            prev = tmp_j->duetime;
        }
    }

    if(prev < problem->H_max) {
        tmp_interval = interval_alloc(prev, problem->H_max, jobarray, njobs);
        g_ptr_array_add(tmp_array, tmp_interval);
    }

    /** calculate the new intervals */
    for(unsigned i = 0; i < tmp_array->len; ++i) {
        tmp_interval = (interval *) g_ptr_array_index(tmp_array, i);
        int count = 0;
        for (size_t j = 0; j < tmp_interval->sigma->len - 1; j++) {
            for (size_t k = j + 1; k < tmp_interval->sigma->len; k++) {
                j1 = (Job *) g_ptr_array_index(tmp_interval->sigma, j);
                j2 = (Job *) g_ptr_array_index(tmp_interval->sigma, k);
                if((double)j1->weight/j1->processingime >= (double)j2->weight/j2->processingime) {
                    if (!check_interval(j1, j2, i, tmp_array)) {
                        count++;
                        check_interval2(j1, j2, i,tmp_array);

                    }
                }
            }
        }
        if (count) {
            printf("count = %d \n", count);
        }
    }

        //
        //intervals_free(tmp_interval);






    CLEAN:
    g_ptr_array_free(tmp_array, TRUE);
    return val;
}
