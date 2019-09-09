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
    int processing_time;
    int release_time;
    int due_time;
    int index;
    int num_layers;
    int *pos_interval;
} Job;

typedef struct _solution {
    PartList *part;
    Job **    perm;
    int *     c;
    int *     u;
    int       tw;
    int       b;
    int       nb_jobs;
    int       nb_machines;
    int       off;
} Solution;

/**
 * Initialization of a solution type
 */
void solution_init(Solution *sol);
/**
 * free all dynamic allocated memory of solution type
 */
void solution_free(Solution **sol);
Solution *solution_alloc(int nb_machines, int nb_jobs, int off);

/**
 * functions to work on solution type data
 */
void solution_print(Solution *sol);
int solution_copy(Solution *dest, Solution *src);
int solution_update(Solution *dest, Solution *src);
int solution_check(PartList *part, int job_count);

Job *job_alloc(int *p, int *w, int *d);
void job_init(Job *job, int p, int d, int w);
void g_job_free(void *set);
void reset_nb_layers(GPtrArray *jobs);
void g_print_job(gpointer data, gpointer user_data);
gint g_job_compare_edd(const void *a, const void *b, void *data);
void g_set_jobarray_job(gpointer data, gpointer user_data);
void g_print_jobarray(gpointer data, gpointer user_data);
void g_print_machine(gpointer data, gpointer user_data);
void g_set_sol_perm(gpointer data, gpointer user_data);
// void g_reset_nb_layers(gpointer data, gpointer user_data);
void g_reset_num_layers(gpointer data, gpointer user_data);

inline int value_Fj(int C, Job *j) { return j->weight * CC_MAX(0, C - j->due_time); }
int value_diff_Fij(int C, Job *i, Job *j);
int bool_diff_Fij(int , Job *, Job *);
int arctime_diff_Fij(int weight, Job* i, Job* j); 
int solution_canonical_order(Solution *sol, GPtrArray *intervals);
void solution_calculate_all(Solution *sol);
void solution_calculate_partition_all(Solution *sol, GPtrArray *intervals);
void solution_calculate_partition_machine(Solution * sol,
                                          GPtrArray *intervals,
                                          int        m);
void solution_calculate_machine(Solution *sol, int m);
int solution_arctime_order(Solution *sol);
#ifdef __cplusplus
}
#endif

#endif  // INCLUDE_DATASTRUCTSOL_H_
