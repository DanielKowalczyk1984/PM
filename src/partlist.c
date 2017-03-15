#include <stdio.h>
#include <datastructsol.h>
#include <util.h>

int partition_order(const void *a, const void *b, void *data) {
    (void)data;
    const int *v1 = (const int *)&(((const Job *)a)->job);
    const int *w1 = (const int *)&(((const Job *)b)->job);
    return (*v1 - *w1);
}

void partlist_free(partlist *part) {
    if (part) {
        if (part->machine != (GPtrArray *)NULL) {
            g_ptr_array_free(part->machine, TRUE);
        }
    }
}

void partlist_init(partlist *part) {
    if (part) {
        part->c = 0;
        part->tw = 0;
        part->machine = g_ptr_array_new();
        part->used = 1;
    }
}

void joblist_init(joblist *jlist) {
    if (jlist) {
        jlist->part = (partlist *)NULL;
    }
}

void partition_init(partlist *part, joblist *jlist, int nbpart, int njobs) {
    int i;

    for (i = 0; i < nbpart; i++) {
        partlist_init(&part[i]);
        part[i].key = i;
    }

    for (i = 0; i < njobs; i++) {
        joblist_init(&jlist[i]);
    }
}
