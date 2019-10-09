#ifndef JOB_H
#define JOB_H

#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include <util.h>

typedef struct _Job {
    int job;
    int weight;
    int processing_time;
    int release_time;
    int due_time;
    int index;
    int num_layers;
    int *pos_interval;
} Job;

inline int value_Fj(int C, Job *j) { return j->weight * CC_MAX(0, C - j->due_time); }
int value_diff_Fij(int C, Job *i, Job *j);
int bool_diff_Fij(int , Job *, Job *);
int arctime_diff_Fij(int weight, Job* i, Job* j);

Job *job_alloc(int *p, int *w, int *d);
void g_job_free(void* set);
void job_init(Job *job, int p, int d, int w);

void g_set_jobarray_job(gpointer data, gpointer user_data);
void g_print_jobarray(gpointer data, gpointer user_data);
void g_print_machine(gpointer data, gpointer user_data);
void g_print_job(gpointer data, gpointer user_data);


#ifdef __cplusplus
}
#endif

#endif // JOB_H
