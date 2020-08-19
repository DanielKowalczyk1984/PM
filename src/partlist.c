#include <partlist.h>
#include <solution.h>
#include <stdio.h>
#include <util.h>

int partition_order(const void* a, const void* b, void* data) {
    (void)data;
    const int* v1 = (const int*)&(((const Job*)a)->job);
    const int* w1 = (const int*)&(((const Job*)b)->job);
    return (*v1 - *w1);
}

void partlist_free(PartList* part) {
    if (part) {
        if (part->machine != (GPtrArray*)NULL) {
            g_ptr_array_free(part->machine, TRUE);
        }
        if (part->Q != (GPtrArray**)NULL) {
            CC_IFFREE(part->Q, GPtrArray*);
        }
        if (part->Q_in != (GPtrArray**)NULL) {
            CC_IFFREE(part->Q_in, GPtrArray*);
        }
    }
}

void partlist_init(PartList* part) {
    if (part) {
        part->c = 0;
        part->tw = 0;
        part->machine = g_ptr_array_new();
        part->used = 1;
    }
}

// void partition_init(PartList* part, int nb_part, MAYBE_UNUSED int nb_jobs) {
//     int i;

//     for (i = 0; i < nb_part; i++) {
//         partlist_init(&part[i]);
//     }
// }
