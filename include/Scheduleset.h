#ifndef SCHEDULESET_H
#define SCHEDULESET_H
#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>

typedef struct Scheduleset {
    int         count;
    int         size;
    int         age;
    int         totweight;
    int         totwct;
    GHashTable *table;
    int *       members;
    int *       C;
    int         id;
} Scheduleset;

void Scheduleset_SWAP(Scheduleset *c1, Scheduleset *c2, Scheduleset *t);

/*Initialization and free memory for the colorset*/
void Scheduleset_init(Scheduleset *set);
void Scheduleset_free(Scheduleset *set);
void Schedulesets_free(Scheduleset **set, int *nsets);
int COLORcopy_sets(Scheduleset **dsts,
                   int *         nsets,
                   Scheduleset * src_s,
                   int           src_nsets);

/*Check if the coloring is feasible*/
int Scheduleset_check_set(Scheduleset *set, int vcount);
int Scheduleset_check(Scheduleset *set, int ccount, int vcount);

/*Transformation of covers*/
int transform_into_coloring(int           vcount,
                            int *         ncolors,
                            Scheduleset **colorclasses);

/*Sorting Schedulesets*/
void Scheduleset_quicksort(Scheduleset *cclasses,
                           int          ccount,
                           int (*functionPtr)(Scheduleset *, Scheduleset *));
void Scheduleset_permquicksort(int *        perm,
                               Scheduleset *cclasses,
                               int          ccount,
                               int (*functionPtr)(Scheduleset *,
                                                  Scheduleset *));
int Scheduleset_less(Scheduleset *c1, Scheduleset *c2);
int Scheduleset_more(Scheduleset *c1, Scheduleset *c2);
int Scheduleset_less_wct(Scheduleset *c1, Scheduleset *c2);
int print_schedule(Scheduleset *cclasses, int ccount);
int Scheduleset_max(Scheduleset *cclasses, int ccount);
int update_Schedulesets(Scheduleset **dst,
                        int *         ndst,
                        Scheduleset * src,
                        int           nsrc);
int add_Schedulesets(Scheduleset **dst, int *ndst, Scheduleset *src, int nsrc);
int Scheduleset_less_totweight(Scheduleset *c1, Scheduleset *c2);
int Scheduleset_more_totweight(Scheduleset *c1, Scheduleset *c2);

#ifdef __cplusplus
}
#endif
#endif  // SCHEDULESET_H
