#ifndef SCHEDULESET_H
#define SCHEDULESET_H
#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include <solution.h>

typedef struct scheduleset {
    int         count;
    int         size;
    int         age;
    int         totweight;
    int         totwct;
    GHashTable *table;
    GPtrArray  *partial_machine;
    int *       members;
    int *       C;
    int         id_right;
    int         id_left;
    int         totwct_right;
    int         totwct_left;
    int         c_right;
    int         c_left;
    int         nb_interval;
} scheduleset;

void scheduleset_SWAP(scheduleset *c1, scheduleset *c2, scheduleset *t);

/*Initialization and free memory for the colorset*/
void scheduleset_init(scheduleset *set);
void scheduleset_free(scheduleset *set);
void schedulesets_free(scheduleset **set, int *nsets);
int COLORcopy_sets(scheduleset **dsts,
                   int *         nsets,
                   scheduleset * src_s,
                   int           src_nsets);

/*Check if the coloring is feasible*/
int scheduleset_check_set(scheduleset *set, int vcount);
int scheduleset_check(scheduleset *set, int ccount, int vcount);

/*Transformation of covers*/
int transform_into_coloring(int           vcount,
                            int *         ncolors,
                            scheduleset **colorclasses);

/*Sorting schedulesets*/
void scheduleset_quicksort(scheduleset *cclasses,
                           int          ccount,
                           int (*functionPtr)(scheduleset *, scheduleset *));
void scheduleset_permquicksort(int *        perm,
                               scheduleset *cclasses,
                               int          ccount,
                               int (*functionPtr)(scheduleset *,
                                                  scheduleset *));
int scheduleset_less(scheduleset *c1, scheduleset *c2);
int scheduleset_more(scheduleset *c1, scheduleset *c2);
int scheduleset_less_wct(scheduleset *c1, scheduleset *c2);
int print_schedule(scheduleset *cclasses, int ccount);
int scheduleset_max(scheduleset *cclasses, int ccount);
int update_schedulesets(scheduleset **dst,
                        int *         ndst,
                        scheduleset * src,
                        int           nsrc);
int add_schedulesets(scheduleset **dst, int *ndst, scheduleset *src, int nsrc);
int scheduleset_less_totweight(scheduleset *c1, scheduleset *c2);
int scheduleset_more_totweight(scheduleset *c1, scheduleset *c2);

/** new approach for columns */
void g_scheduleset_free(void *set);
scheduleset *scheduleset_create(Job **job_array, int nbjobs, GPtrArray *intervals, int u);
scheduleset *scheduleset_from_solution(Job **job_array, int begin, int end, GPtrArray *intervals, int u);

#ifdef __cplusplus
}
#endif
#endif  // SCHEDULESET_H
