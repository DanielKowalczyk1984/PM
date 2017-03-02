#include <stdio.h>
#include "datastructsol.h"
#include "util.h"

int partition_order(const void *a, const void *b, void *data) {
    (void) data;
    const int *v1 = (const int *) & (((const Job *)a)->job);
    const int *w1 = (const int *) & (((const Job *)b)->job);
    return (*v1 - *w1);
}

void partlist_free(partlist *part) {
    if (part) {
        if (part->list != (GQueue *) NULL) {
            g_queue_free(part->list);
        }
        if(part->machine != (GPtrArray *) NULL) {
            g_ptr_array_free(part->machine, TRUE);
        }

        part->list = (GQueue *) NULL;
    }
}

void partlist_init(partlist *part) {
    if (part) {
        part->list = g_queue_new();
        part->c = 0;
        part->tw = 0;
        part->machine = g_ptr_array_new();
        part->used = 1;
    }
}

void joblist_init(joblist *jlist) {
    if (jlist) {
        jlist->part = (partlist *) NULL;
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


int partlist_insert(partlist *part, joblist *jlist, Job *job) {
    int val = 0;

    if (jlist[job->job].part != NULL) {
        fprintf(stderr, "Error: double insertion\n");
        val = 1;
        goto CLEAN;
    }

    jlist[job->job].part = part;
    g_queue_push_tail(part->list, job);
    part->c += job->processingime;
    part->tw += job->weight*CC_MAX(part->c - job->duetime, 0);
CLEAN:
    return val;
}


int partlist_delete(joblist *jlist, Job *job) {
    int val = 0;
    partlist *p = NULL;

    if (jlist[job->job].part == NULL) {
        fprintf(stderr, "Error deleting a job that is not assigned\n");
        val = 1;
        goto CLEAN;
    }

    p = jlist[job->job].part;
    g_queue_remove(p->list, job);
    p->c -= job->processingime;
    jlist[job->job].part = NULL;
CLEAN:
    return val;
}

void partlist_move(partlist *part, joblist *jlist, Job *job) {
    if (jlist[job->job].part != NULL) {
        partlist_delete(jlist, job);
        partlist_insert(part, jlist, job);
    } else {
        partlist_insert(part, jlist, job);
    }
}
