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
} partlist;

typedef struct _joblist { partlist *part; } joblist;

void partlist_free(partlist *part);
void partlist_init(partlist *part);
void joblist_init(joblist *vlist);
void partition_init(partlist *part, joblist *vlist, int nbpart, int jcount);
int partition_order(const void *a, const void *b, void *data);
int find_vertex(const void *a, const void *b);

void partlist_permquicksort(int *     perm,
                            partlist *part,
                            int       nbpart,
                            int (*functionPtr)(partlist *, partlist *));
int partlist_more_totweight(partlist *c1, partlist *c2);
#ifdef __cplusplus
}
#endif
#endif  // PARTLIST_H
