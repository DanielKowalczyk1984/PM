#include <wct.h>

static int _job_compare_edd(const void *a, const void *b);

static int _job_compare_edd(const void *a, const void *b) {
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
    int      njobs = problem->njobs;
    Job **   _ojobarray = (Job **)NULL;
    wctdata *pd = (wctdata *)NULL;

    /** Initialize jobarray of rootnode */
    for (int i = 0; i < njobs; ++i) {
        problem->psum += problem->jobarray[i].processingime;
        problem->pmax =
            CC_MAX(problem->pmax, problem->jobarray[i].processingime);
        problem->pmin =
            CC_MIN(problem->pmin, problem->jobarray[i].processingime);
        problem->dmax = CC_MAX(problem->dmax, problem->jobarray[i].duetime);
        problem->dmin = CC_MIN(problem->pmin, problem->jobarray[i].duetime);
    }

    pd = &(problem->root_pd);
    /** Calculate H_max */
    temp = problem->psum - problem->pmax;
    temp_dbl = (double)temp;
    temp_dbl = temp_dbl / problem->nmachines + problem->pmax;
    problem->T = pd->H_max = (int)ceil(temp_dbl);
    printf("H_max = %d\n", pd->H_max);
    /** Create edd ordered array*/
    _ojobarray = CC_SAFE_MALLOC(njobs, Job *);
    CCcheck_NULL_2(_ojobarray, "Failed to allocate memory");

    for (int i = 0; i < njobs; ++i) {
        _ojobarray[i] = &(problem->jobarray[i]);
    }

    qsort((void *)_ojobarray, njobs, sizeof(Job *), _job_compare_edd);

    for (int i = 0; i < problem->njobs; ++i) {
        _ojobarray[i]->job = i;
    }

    pd->njobs = problem->njobs;
    pd->jobarray = problem->jobarray;
    pd->nmachines = problem->nmachines;
    problem->ojobarray = _ojobarray;
CLEAN:

    if (val) {
        CC_IFFREE(_ojobarray, Job *);
    }

    return val;
}
