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
    int         totweight;
    int         totwct;
    GHashTable *table;
    GPtrArray  *jobs;
    GPtrArray *e_list;
    int        *nb;
    int         id;
} scheduleset;

void scheduleset_SWAP(scheduleset *c1, scheduleset *c2, scheduleset *t);
scheduleset *scheduleset_alloc_bis(int nbjobs) ;
void scheduleset_init_bis(scheduleset *set);

/*Initialization and free memory for the colorset*/
void scheduleset_init(scheduleset *set);
void scheduleset_free(scheduleset *set);
void schedulesets_free(scheduleset **set, int *nsets);

/*Sorting schedulesets*/
int scheduleset_less(scheduleset *c1, scheduleset *c2);
int print_schedule(scheduleset *cclasses, int ccount);
int scheduleset_max(scheduleset *cclasses, int ccount);
gint g_scheduleset_less(gconstpointer a,gconstpointer b);
void g_scheduleset_print(gpointer data, gpointer user_data);

/** new approach for columns */
void g_scheduleset_free(void *set);
void g_sum_processing_time(gpointer data, gpointer user_data);
scheduleset *scheduleset_from_solution(GPtrArray *machine, int nbjobs);
scheduleset *scheduleset_alloc(int nbjobs);

#ifdef __cplusplus
}
#endif
#endif  // SCHEDULESET_H
