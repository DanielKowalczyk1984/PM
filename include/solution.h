#ifndef INCLUDE_DATASTRUCTSOL_H_
#define INCLUDE_DATASTRUCTSOL_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include <partlist.h>

typedef struct _Job {
    int job;
    int weight;
    int processingime;
    int releasetime;
    int duetime;
    int index;
} Job;

typedef struct _solution {
    partlist *part;
    Job **    perm;
    int *     c;
    int       tw;
    int       b;
    int       njobs;
    int       nmachines;
    int       off;
} solution;

void solution_init(solution *sol);
void solution_free(solution **sol);
solution *solution_alloc(int nmachines, int njobs, int off);

void solution_print(solution *sol);
int solution_copy(solution *dest, solution *src);
int solution_update(solution *dest, solution *src);
int solution_check(partlist *part, int jcount);

#ifdef __cplusplus
}
#endif

#endif  // INCLUDE_DATASTRUCTSOL_H_
