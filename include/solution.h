#ifndef INCLUDE_DATASTRUCTSOL_H_
#define INCLUDE_DATASTRUCTSOL_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include <partlist.h>
#include <util.h>

typedef struct _Job {
    int job;
    int weight;
    int processingime;
    int releasetime;
    int duetime;
    int index;
    int nb_layers;
    int *pos_interval;
} Job;

typedef struct _solution {
    partlist *part;
    Job **    perm;
    int *     c;
    int *     u;
    int       tw;
    int       b;
    int       njobs;
    int       nmachines;
    int       off;
} solution;

/**
 * Initialization of a solution type
 */
void solution_init(solution *sol);
/**
 * free all dynamic allocated memory of solution type
 */
void solution_free(solution **sol);
solution *solution_alloc(int nmachines, int njobs, int off);

/**
 * functions to work on solution type data
 */
void solution_print(solution *sol);
int solution_copy(solution *dest, solution *src);
int solution_update(solution *dest, solution *src);
int solution_check(partlist *part, int jcount);

Job *job_alloc(int *p, int *w, int *d);
void g_job_free(void *set);
void reset_nblayers(GPtrArray *jobs);
void g_print_job(gpointer data, gpointer user_data);
gint g_job_compare_edd(const void *a, const void *b, void *data);
void g_set_jobarray_job(gpointer data, gpointer user_data);
void g_print_jobarray(gpointer data, gpointer user_data);
void g_print_machine(gpointer data, gpointer user_data);
void g_set_sol_perm(gpointer data, gpointer user_data);
void g_reset_nb_layers(gpointer data, gpointer user_data);

inline int value_Fj(int C, Job *j) { return j->weight * CC_MAX(0, C - j->duetime); }
int value_diff_Fij(int C, Job *i, Job *j);

int solution_canonical_order(solution *sol, GPtrArray *intervals);
void solution_calculate_all(solution *sol);
void solution_calculate_partition_all(solution *sol, GPtrArray *intervals);
void solution_calculate_partition_machine(solution * sol,
                                          GPtrArray *intervals,
                                          int        m);
void solution_calculate_machine(solution *sol, int m);

#ifdef __cplusplus
}
#endif

#endif  // INCLUDE_DATASTRUCTSOL_H_
