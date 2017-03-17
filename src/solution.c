#include <solution.h>
#include <string.h>
#include <util.h>

gint comparefunc(const void *a, const void *b, void *data);
gint compare_func(const void *a, const void *b);
gint order_weight(gconstpointer a, gconstpointer b, void *data);

void solution_init(solution *sol) {
    if (sol) {
        sol->part = (partlist *)NULL;
        sol->perm = (Job **)NULL;
        sol->c = (int *)NULL;
        sol->nmachines = 0;
        sol->njobs = 0;
        sol->tw = 0;
        sol->b = 0;
        sol->off = 0;
    }
}

void solution_free(solution **sol) {
    if (*sol) {
        for (int i = 0; i < (*sol)->nmachines; ++i) {
            partlist_free((*sol)->part + i);
        }

        CC_IFFREE((*sol)->part, partlist);
        CC_IFFREE((*sol)->perm, Job *);
        CC_IFFREE((*sol)->c, int);
        CC_IFFREE((*sol), solution);
    }
}

solution *solution_alloc(int nmachines, int njobs, int off) {
    int       val = 0;
    int       i;
    solution *sol = CC_SAFE_MALLOC(1, solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory");
    solution_init(sol);
    sol->nmachines = nmachines;
    sol->njobs = njobs;
    sol->tw = 0;
    sol->b = 0;
    sol->off = off;
    sol->part = CC_SAFE_MALLOC(nmachines, partlist);
    CCcheck_NULL_2(sol->part, "Failed to allocate memory to part");

    for (i = 0; i < nmachines; ++i) {
        partlist_init(sol->part + i);
        (sol->part + i)->key = i;
    }

    sol->perm = CC_SAFE_MALLOC(njobs, Job *);
    CCcheck_NULL_2(sol->perm, "Failed to allocate memory to perm");
    sol->c = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(sol->c, "Failed to allocate memory");
    fill_int(sol->c, sol->njobs, 0);

    for (i = 0; i < njobs; ++i) {
        sol->perm[i] = (Job *)NULL;
    }

CLEAN:

    if (val) {
        solution_free(&sol);
    }

    return sol;
}

gint comparefunc(const void *a, const void *b, void *data) {
    (void)data;
    const int *v = &(((const Job *)a)->job);
    const int *w = &(((const Job *)b)->job);
    return *v - *w;
}

gint order_weight(gconstpointer a, gconstpointer b, void *data) {
    (void)data;
    const int *v = &(((const Job *)a)->weight);
    const int *w = &(((const Job *)b)->weight);
    return -(*v - *w);
}

static void print_machine(gpointer j, gpointer data) {
    Job *tmp = (Job *)j;
    printf("%3d ", tmp->job);
}

void solution_print(solution *sol) {
    for (int i = 0; i < sol->nmachines; ++i) {
        printf("Machine %-1d: ", sol->part[i].key);
        g_ptr_array_foreach(sol->part[i].machine, print_machine, NULL);
        printf("with C =  %d, wC = %d and %u jobs\n", sol->part[i].c,
               sol->part[i].tw, sol->part[i].machine->len);
    }

    printf("with total weighted tardiness %d\n", sol->tw + sol->off);
}

int solution_copy(solution *dest, solution *src) {
    int val = 0;
    dest = solution_alloc(src->nmachines, src->njobs, src->off);
    CCcheck_val_2(val, "Failed in  solution_alloc");
    dest->tw = src->tw;
    dest->b = src->b;
    dest->off = src->off;

    for (int i = 0; i < dest->nmachines; i++) {
        dest->part[i].key = src->part[i].key;
        dest->part[i].tw = src->part[i].tw;
        dest->part[i].c = src->part[i].c;
    }

CLEAN:

    if (val) {
        solution_free(&dest);
    }

    return val;
}

int solution_update(solution *dest, solution *src) {
    int val = 0;
    dest->tw = src->tw;
    dest->b = src->b;
    dest->nmachines = src->nmachines;
    dest->njobs = src->njobs;
    dest->off = src->off;

    for (int i = 0; i < dest->nmachines; i++) {
        g_ptr_array_free(dest->part[i].machine, TRUE);
        dest->part[i].machine = g_ptr_array_new();

        for (unsigned j = 0; j < src->part[i].machine->len; ++j) {
            g_ptr_array_add(dest->part[i].machine,
                            g_ptr_array_index(src->part[i].machine, j));
        }

        dest->part[i].tw = src->part[i].tw;
        dest->part[i].c = src->part[i].c;
    }

    memcpy(dest->perm, src->perm, src->njobs * sizeof(Job *));
    memcpy(dest->c, src->c, dest->njobs * sizeof(int));
    return val;
}

void partlist_permquicksort(int *     perm,
                            partlist *part,
                            int       nbpart,
                            int (*functionPtr)(partlist *, partlist *)) {
    int      i, j, temp;
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
