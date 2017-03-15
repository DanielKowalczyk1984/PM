////////////////////////////////////////////////////////////////
//                                                            //
//  allocr.c                                                  //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <util.h>

void *CCutil_allocrus(size_t size) {
    void *mem = (void *)NULL;

    if (size == 0) {
        fprintf(stderr, "Warning: 0 bytes allocated\n");
        return mem;
    }

    mem = (void *)malloc(size);

    if (mem == (void *)NULL) {
        fprintf(stderr, "Out of memory. Asked for %d bytes\n at %d in file %s",
                (int)size, __LINE__, __FILE__);
    }

    return mem;
}

void CCutil_freerus(void *ptr) {
    if (!ptr) {
        fprintf(stderr, "Warning: null pointer freed\n");
        return;
    }

    free(ptr);
}

void *CCutil_reallocrus(void *ptr, size_t size) {
    if (!ptr) {
        return CCutil_allocrus(size);
    } else {
        void *newptr;
        newptr = (void *)realloc(ptr, size);

        if (!newptr) {
            fprintf(stderr, "Out of memory.  Tried to grow to %d bytes\n",
                    (int)size);
        }

        return newptr;
    }
}

int CCutil_reallocrus_scale(
    void **pptr, int *pnnum, int count, double scale, size_t size) {
    int   newsize = (int)(((double)*pnnum) * scale);
    void *ptr;

    if (newsize < *pnnum + 1000) {
        newsize = *pnnum + 1000;
    }

    if (newsize < count) {
        newsize = count;
    }

    ptr = CCutil_reallocrus(*pptr, newsize * size);

    if (!ptr) {
        return 1;
    } else {
        *pptr = ptr;
        *pnnum = newsize;
        return 0;
    }
}

int CCutil_reallocrus_count(void **pptr, int count, size_t size) {
    void *ptr = CCutil_reallocrus(*pptr, count * size);

    if (!ptr) {
        return 1;
    } else {
        *pptr = ptr;
        return 0;
    }
}
