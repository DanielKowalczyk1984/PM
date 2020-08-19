#include <job.h>

void g_print_job(gpointer data, MAYBE_UNUSED gpointer user_data) {
    Job* a = (Job*)data;
    printf("%d ", a->job);
}

void g_job_free(void* set) {
    Job* tmp = (Job*)set;
    if (tmp) {
        CC_IFFREE(tmp->pos_interval, int);
        CC_IFFREE(tmp, Job);
    }
}

Job* job_alloc(int* p, int* w, int* d) {
    Job* j = CC_SAFE_MALLOC(1, Job);
    j->processing_time = *p;
    j->due_time = *d;
    j->weight = *w;
    j->num_layers = 0;
    j->pos_interval = (int*)NULL;
    return j;
}

void job_init(Job* job, int p, int w, int d) {
    job->processing_time = p;
    job->weight = w;
    job->due_time = d;
}

void g_set_jobarray_job(gpointer data, gpointer user_data) {
    Job* j = (Job*)data;
    int* i = (int*)user_data;
    j->job = *i;
    j->num_layers = 0;
    (*i)++;
}

void g_print_jobarray(gpointer data, MAYBE_UNUSED gpointer user_data) {
    Job* j = (Job*)data;
    g_print("Job %d: %d %d %d %f\n", j->job, j->processing_time, j->due_time,
            j->weight, (double)j->weight / j->processing_time);
}

void g_print_machine(gpointer data, MAYBE_UNUSED gpointer user_data) {
    Job* j = (Job*)data;
    g_print("%d ", j->job);
}

extern inline int value_Fj(int C, Job* j);

int value_diff_Fij(int C, Job* i, Job* j) {
    int val = value_Fj(C + i->processing_time - j->processing_time, i);
    val += value_Fj(C + i->processing_time, j);
    val -= value_Fj(C, j);
    val -= value_Fj(C + i->processing_time, i);
    return val;
}

int bool_diff_Fij(int weight, Job* _prev, Job* tmp_j) {
    return (_prev == NULL) ? 1
                           : (value_diff_Fij(weight + tmp_j->processing_time,
                                             _prev, tmp_j) >= 0);
}

int arctime_diff_Fij(int weight, Job* i, Job* j) {
    int val = value_Fj(weight, i);
    val += value_Fj(weight + j->processing_time, j);
    val -= value_Fj(weight - i->processing_time + j->processing_time, j);
    val -= value_Fj(weight + j->processing_time, i);
    return val;
}
