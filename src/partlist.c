#include "partlist.h"
#include "util.h"

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
