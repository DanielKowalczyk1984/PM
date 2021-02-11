#include <interval.h>
#include <wct.h>
#include "job.h"

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
