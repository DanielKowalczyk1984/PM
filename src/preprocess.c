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
    int i = 0;
    int      temp = 0;
    double   temp_dbl = 0.0;
    wctdata *pd = (wctdata *)NULL;

    g_ptr_array_foreach(problem->g_job_array, g_problem_summary_init, problem);

    pd = &(problem->root_pd);
    pd->njobs = problem->njobs;
    pd->nmachines = problem->nmachines;
    /** Calculate H_max */
    temp = problem->psum - problem->pmax;
    temp_dbl = (double)temp;
    temp_dbl = temp_dbl / problem->nmachines + problem->pmax;
    problem->H_max = pd->H_max = (int)ceil(temp_dbl);
    printf("H_max = %d\n", problem->H_max);

    g_ptr_array_sort_with_data(problem->g_job_array, g_job_compare_edd, NULL);
    g_ptr_array_foreach(problem->g_job_array, g_set_jobarray_job, &i);

    find_division(problem);
CLEAN:

    return val;
}



int find_division(wctproblem *problem){
    int val = 0;
    int tmp;
    interval *tmp_interval_ptr;
    int prev;
    int njobs = problem->njobs;
    GPtrArray *tmp_e = g_ptr_array_new_with_free_func(intervals_free);
    GPtrArray* jobarray = problem->g_job_array;

    /** Find initial partition */
    prev = 0;
    for(unsigned i = 0; i < njobs && prev != problem->H_max; ++i) {
        tmp = CC_MIN(problem->H_max, ((Job *)jobarray->pdata[i])->duetime);
        if(prev < tmp) {
            tmp_interval_ptr = interval_alloc(prev, ((Job *)jobarray->pdata[i])->duetime, jobarray, njobs);
            CCcheck_NULL_2(tmp_interval_ptr, "Failed to allocate memory")
            prev = ((Job *)jobarray->pdata[i])->duetime;
            g_ptr_array_add(tmp_e, tmp_interval_ptr);
        }
    }

    if(prev < problem->H_max) {
        tmp_interval_ptr = interval_alloc(prev, problem->H_max, jobarray, njobs);
        g_ptr_array_add(tmp_e, tmp_interval_ptr);
    }

    g_ptr_array_foreach(tmp_e, g_print_interval, NULL);



    CLEAN:
    g_ptr_array_free(tmp_e, TRUE);
    return val;
}
