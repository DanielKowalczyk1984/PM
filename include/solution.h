#ifndef INCLUDE_DATASTRUCTSOL_H_
#define INCLUDE_DATASTRUCTSOL_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include <partlist.h>

typedef struct _Job {
    int job;
    int weight;
    int processingime;
    int releasetime;
    int duetime;
    int index;
} Job;

typedef struct _solution {
    partlist *part;
    Job **    perm;
    int *     c;
    int       tw;
    int       b;
    int       njobs;
    int       nmachines;
    int       off;
} solution;

void solution_init(solution *sol);
void solution_free(solution **sol);
solution *solution_alloc(int nmachines, int njobs, int off);

void solution_print(solution *sol);
int solution_copy(solution *dest, solution *src);
int solution_update(solution *dest, solution *src);
int solution_check(partlist *part, int jcount);

Job *job_alloc(int *p, int *w, int *d);
void g_print_job(gpointer data, gpointer user_data);
gint g_job_compare_edd(const void *a, const void *b, void *data);
void g_set_jobarray_job(gpointer data, gpointer user_data);
void g_print_jobarray(gpointer data, gpointer user_data);
void g_set_sol_perm(gpointer data, gpointer user_data);

#ifdef __cplusplus
}
#endif

#endif  // INCLUDE_DATASTRUCTSOL_H_
