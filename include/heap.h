#ifndef _HEAP_H
#define _HEAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#define CCutil_MAXDBL (1e30)
#define CCutil_MAXINT (2147483647)
typedef struct heapelm {
    int   key;
    void* obj;
} HeapElement;

typedef struct Heap_t {
    int  end;
    int  size;
    int* perm;
    int* iperm;

    HeapElement* elms;
} HeapContainer;

int heapcontainer_init(HeapContainer**heap, int size),
    heapcontainer_free(HeapContainer*heap),
    heapcontainer_free_all(HeapContainer*heap),
    heapcontainer_insert(HeapContainer*heap, int key, void*obj),
    heapcontainer_remove(HeapContainer*heap, int href),
    heapcontainer_get_key(const HeapContainer*heap, int href),
    heapcontainer_size(const HeapContainer*heap),
    heapcontainer_decrease_key(HeapContainer*heap, int href, int new_key),
    heapcontainer_relabel(HeapContainer*heap, int href, int new_key);
void* heapcontainer_get_obj(const HeapContainer* heap, int href);

void *heapcontainer_min(HeapContainer*heap),
    heapcontainer_reset(HeapContainer*heap),
    heapcontainer_reset_free(HeapContainer*heap);

#ifdef __cplusplus
}
#endif
#endif
