#ifndef PARTLIST_H
#define PARTLIST_H

#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>

typedef struct _partlist {
    GPtrArray*  machine;
    int         c;
    int         tw;
    int         used;
    GPtrArray** Q_in;
    GPtrArray** Q;
} PartList;

void partlist_free(PartList* part);
void partlist_init(PartList* part);
#ifdef __cplusplus
}
#endif
#endif  // PARTLIST_H
