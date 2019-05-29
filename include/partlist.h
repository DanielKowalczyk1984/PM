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
    int        key;
    int        used;
} PartList;

typedef struct _joblist { PartList *part; } joblist;

void partlist_free(PartList *part);
void partlist_init(PartList *part);
void joblist_init(joblist *vlist);
void partition_init(PartList *part, joblist *vlist, int nbpart, int jcount);
int partition_order(const void *a, const void *b, void *data);
int find_vertex(const void *a, const void *b);

void partlist_permquicksort(int *     perm,
                            PartList *part,
                            int       nbpart,
                            int (*functionPtr)(PartList *, PartList *));
int partlist_more_totweight(PartList *c1, PartList *c2);
#ifdef __cplusplus
}
#endif
#endif  // PARTLIST_H
