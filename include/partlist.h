#ifndef PARTLIST_H
#define PARTLIST_H

#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>

typedef struct _partlist {
    GPtrArray *machine;
    int        c;
    int        tw;
    int        used;
    GPtrArray **Q_in;
    GPtrArray **Q;
} PartList;

void partlist_free(PartList *part);
void partlist_init(PartList *part);
// void joblist_init(joblist *vlist);
void partition_init(PartList *part, int nb_part, int nb_jobs);
int partition_order(const void *a, const void *b, void *data);
int find_vertex(const void *a, const void *b);

void partlist_permquicksort(int *     perm,
                            PartList *part,
                            int       nb_part,
                            int (*functionPtr)(PartList *, PartList *));
int partlist_more_totweight(PartList *c1, PartList *c2);
#ifdef __cplusplus
}
#endif
#endif  // PARTLIST_H
