#include "datastructsol.h"
#include "util.h"
#include <string.h>

gint comparefunc(const void *a, const void *b, void *data);
gint compare_func(const void *a, const void *b);
gint order_weight(gconstpointer a, gconstpointer b, void *data);


void solution_init(solution *sol) {
    if (sol) {
        sol->part     = (partlist *)NULL;
        sol->vlist    = (joblist *)NULL;
        sol->perm     = (Job **)NULL;
        sol->c        = (int*) NULL;
        sol->nmachines   = 0;
        sol->njobs   = 0;
        sol->tw = 0;
        sol->b = 0;
    }
}

void solution_free(solution *sol) {
    int i;

    if (sol) {
        for (i = 0; i < sol->nmachines; ++i) {
            partlist_free(sol->part + i);
        }

        CC_IFFREE(sol->part, partlist);
        CC_IFFREE(sol->vlist, joblist);
        CC_IFFREE(sol->perm, Job*);
        CC_IFFREE(sol->c, int);
        sol->nmachines   = 0;
        sol->tw = 0;
        sol->njobs   = 0;
        sol->b = 0;
    }
}

solution* solution_alloc(int nmachines, int njobs) {
    int val = 0;
    int i;
    
    solution *sol = (solution *) NULL;
    sol = CC_SAFE_MALLOC(1, solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory");
    solution_init(sol);

    sol->nmachines  = nmachines;
    sol->njobs = njobs;
    sol->tw = 0;
    sol->b = 0;

    sol->part = CC_SAFE_MALLOC(nmachines, partlist);
    CCcheck_NULL_2(sol->part, "Failed to allocate memory to part");

    for (i = 0; i < nmachines; ++i) {
        partlist_init(sol->part + i);
        (sol->part + i)->key = i;
    }

    sol->vlist = CC_SAFE_MALLOC(njobs, joblist);
    CCcheck_NULL_2(sol->vlist, "Failed to allocate memory to vlist");
    sol->perm = CC_SAFE_MALLOC(njobs, Job*);
    CCcheck_NULL_2(sol->perm, "Failed to allocate memory to perm");
    sol->c = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(sol->c, "Failed to allocate memory");
    fill_int(sol->c, sol->njobs, 0);

    for (i = 0; i < njobs; ++i) {
        joblist_init(sol->vlist + i);
        sol->perm[i] = (Job *) NULL;
    }

CLEAN:

    if (val) {
        solution_free(sol);
    }

    return sol;
}

gint comparefunc(const void *a, const void *b, void *data) {
    (void) data;
    const int *v = &(((const Job *)a)->job);
    const int *w = &(((const Job *)b)->job);
    return *v - *w;
}

gint compare_func(const void *a, const void *b) {
    const int *v = &(((const partlist *)a)->c);
    const int *w = &(((const partlist *)b)->c);

    if (*v != *w) {
        return *v - *w;
    } else {
        if (*v == 0 || *w == 0) {
            return *v - *w;
        }

        const int *vv = &(((Job *)((const partlist *)a)->list->head->data)->job);
        const int *ww = &(((Job *)((const partlist *)b)->list->head->data)->job);
        return *vv - *ww;
    }
}

gint order_weight(gconstpointer a, gconstpointer b, void *data) {
    (void) data;
    const int *v = &(((const Job *)a)->weight);
    const int *w = &(((const Job *)b)->weight);
    return -(*v - *w);
}

void solution_max(solution *sol) {
    int max = 0;

    for (int i = 0; i < sol->nmachines; ++i) {
        if (sol->part[i].c > max) {
            max = sol->part[i].c;
        }
    }

    sol->tw = max;
}

void solution_unique(solution *sol) {
    int i;
    int nmachines = sol->nmachines;
    GList *it = (GList *) NULL;
    partlist *temp_partlist = (partlist *) NULL;
    int counter = 0;
    qsort(sol->part, sol->nmachines, sizeof(partlist), compare_func);

    /** Compute permutation */
    if (sol->perm == NULL) {
        sol->perm = CC_SAFE_MALLOC(sol->njobs, Job*);
    }

    for (i = 0; i < nmachines; i++) {
        temp_partlist = sol->part + i;
        temp_partlist->key = i;

        for (it = temp_partlist->list->head; it; it = it->next) {
            sol->perm[counter] = ((Job *)it->data);
            sol->vlist[((Job *)it->data)->job].part = sol->part + i;
            counter++;
        }
    }
}

void solution_print(solution *sol) {
    for (int i = 0; i < sol->nmachines; ++i) {
        printf("Machine %d: ", sol->part[i].key);

        for (GList *it = sol->part[i].list->head; it; it = it->next) {
            printf("%d ", ((Job *)it->data)->job);
        }

        printf("with C =  %d, wC = %d and %d jobs\n", sol->part[i].c,
               sol->part[i].tw
               , g_queue_get_length(sol->part[i].list));
    }

    printf("with total weighted tardiness %d\n", sol->tw);
}

void test_SOLUTION(solution *sol) {
    int i;

    for (i = 0; i <  sol->njobs; i++) {
        printf("%d ", i);

        for (GList *it = sol->vlist[i].part->list->head; it; it = it->next) {
            printf("%d ", ((Job *)it->data)->job);
        }

        printf("\n");
    }
}

int solution_copy(solution *dest, solution *src) {
    int val = 0;
    int counter = 0;
    dest = solution_alloc(src->nmachines, src->njobs);
    CCcheck_val_2(val, "Failed in  solution_alloc");
    dest->tw = src->tw;
    dest->b = src->b;

    for (int i = 0; i < dest->nmachines; i++) {
        dest->part[i].key = src->part[i].key;
        dest->part[i].tw = src->part[i].tw;
        g_queue_free(dest->part[i].list);
        dest->part[i].list = (GQueue *) NULL;
        dest->part[i].list = g_queue_copy(src->part[i].list);
        dest->part[i].c = src->part[i].c;

        for (GList *it = dest->part[i].list->head; it; it = it->next) {
            dest->perm[counter] = ((Job *)it->data);
            dest->vlist[((Job *)it->data)->job].part = dest->part + i;
            counter++;
        }
    }

CLEAN:

    if (val) {
        solution_free(dest);
        CC_IFFREE(dest, solution);
    }

    return val;
}

int solution_update(solution *dest, solution *src) {
    int val = 0;
    dest->tw = src->tw;
    dest->b = src->b;
    dest->nmachines   = src->nmachines;
    dest->njobs   = src->njobs;

    for (int i = 0; i < dest->nmachines; i++) {
        g_queue_free(dest->part[i].list);
        dest->part[i].tw = src->part[i].tw;
        dest->part[i].list = g_queue_copy(src->part[i].list);
        dest->part[i].c = src->part[i].c;
    }

    memcpy(dest->perm, src->perm, src->njobs*sizeof(Job*));
    memcpy(dest->vlist, src->vlist,src->njobs*sizeof(joblist));

    return val;
}

void solution_calc(solution *sol, Job *jobarray) {
    int i;
    int nmachines = sol->nmachines;
    GList *it = (GList *) NULL;
    Job *temp_job = (Job *) NULL;
    partlist *temp_partlist = (partlist *) NULL;
    sol->tw = 0;
    sol->b = 0;

    /** Order in WSPT order and compute objective value of this solution */
    for (i = 0; i < nmachines; ++i) {
        temp_partlist = sol->part + i;
        temp_partlist->c = 0;
        g_queue_sort(temp_partlist->list, (GCompareDataFunc)comparefunc, NULL);

        for (it = temp_partlist->list->head; it; it = g_list_next(it)) {
            temp_job = ((Job *)it->data);
            temp_partlist->c += temp_job->processingime ;
            temp_partlist->tw += temp_job->weight *
                                            temp_partlist->c;
            sol->tw += temp_partlist->c * temp_job->weight;
        }
    }

    /** Compute sum weights and sum times */


}

void solution_wct(solution *sol) {
    int i;
    int nmachines = sol->nmachines;
    GList *it = (GList *) NULL;
    Job *temp_job = (Job *) NULL;
    partlist *temp_partlist = (partlist *) NULL;
    sol->tw = 0;
    sol->b = 0;

    /** Order in WSPT order and compute objective value of this solution */
    for (i = 0; i < nmachines; ++i) {
        temp_partlist = sol->part + i;
        temp_partlist->c = 0;
        g_queue_sort(temp_partlist->list, (GCompareDataFunc)comparefunc, NULL);

        for (it = temp_partlist->list->head; it; it = g_list_next(it)) {
            temp_job = ((Job *)it->data);
            temp_partlist->c += temp_job->processingime;
            sol->tw += temp_partlist->c * temp_job->weight;
        }
    }
}

void partlist_permquicksort(int *perm, partlist *part, int nbpart,
                            int (*functionPtr)(partlist *, partlist *)) {
    int i, j, temp;
    partlist t;

    if (nbpart <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(nbpart - 1) / 2], temp);
    i = 0;
    j = nbpart;
    memcpy(&t, &(part[perm[0]]), sizeof(partlist));

    while (1) {
        do {
            i++;
        } while (i < nbpart && (*functionPtr)(&(part[perm[i]]), &t));

        do {
            j--;
        } while ((*functionPtr)(&t, &(part[perm[j]])));

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    partlist_permquicksort(perm, part, j, (*functionPtr));
    partlist_permquicksort(perm + i, part, nbpart - i, (*functionPtr));
}



