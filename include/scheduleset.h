#ifndef SCHEDULESET_H
#define SCHEDULESET_H
#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include <solution.h>
#include <interval.h>

typedef struct scheduleset {
    int         age;
    int         total_processing_time;
    int         total_weighted_completion_time;
    GHashTable* table;
    GPtrArray*  job_list;
    GPtrArray*  edge_list;
    int*        num;
    int         id;
} ScheduleSet;

void scheduleset_SWAP(ScheduleSet *c1, ScheduleSet *c2, ScheduleSet *t);
ScheduleSet *scheduleset_alloc_bis(int nbjobs) ;
void scheduleset_init_bis(ScheduleSet *set);

/*Initialization and free memory for the colorset*/
void scheduleset_init(ScheduleSet *set);
void scheduleset_free(ScheduleSet *set);
void schedulesets_free(ScheduleSet **set, int *nsets);

/*Sorting schedulesets*/
int scheduleset_less(ScheduleSet *c1, ScheduleSet *c2);
int print_schedule(ScheduleSet *cclasses, int ccount);
int scheduleset_max(ScheduleSet *cclasses, int ccount);
gint g_scheduleset_less(gconstpointer a,gconstpointer b);
void g_scheduleset_print(gpointer data, gpointer user_data);

/** new approach for columns */
void g_scheduleset_free(void *set);
void g_sum_processing_time(gpointer data, gpointer user_data);
void g_compute_nblayers_schedule(gpointer data, gpointer user_data);
ScheduleSet *scheduleset_from_solution(GPtrArray *machine, int nbjobs);
ScheduleSet *scheduleset_alloc(int nbjobs);

#ifdef __cplusplus
}
#endif
#endif  // SCHEDULESET_H
