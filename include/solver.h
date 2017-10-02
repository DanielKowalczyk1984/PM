#ifndef SOLVER_H
#define SOLVER_H

#include <wctparms.h>
#include <solution.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PricerSolver PricerSolver;
PricerSolver *newSolver(GPtrArray *interval_list, int njobs, int **sum_p);
PricerSolver *newSolverDP(GPtrArray *interval_list, int njobs, int **sum_p);
PricerSolver *copySolver(PricerSolver *src);
int init_tables(PricerSolver *solver);
int calculate_table(PricerSolver *solver, wctparms *parms);

void copy_solver(PricerSolver **dest, PricerSolver *src);

void freeSolver(PricerSolver *src);
void deletePricerSolver(PricerSolver *solver);

int add_conflict_constraints(PricerSolver *solver,
                             wctparms *    parms,
                             int *         elist_same,
                             int           ecount_same,
                             int *         elist_differ,
                             int           ecount_differ);
int free_conflict_constraints(PricerSolver *solver,
                              wctparms *    parms,
                              int           ecount_same,
                              int           ecount_differ);
int add_one_conflict(
    PricerSolver *solver, wctparms *parms, Job *v1, Job *v2, int same);

void iterate_dd(PricerSolver *solver);
void iterate_zdd(PricerSolver *solver);

void set_release_due_time(PricerSolver *solver, Job *jobarray);

size_t get_datasize(PricerSolver *solver);
size_t get_numberrows_zdd(PricerSolver *solver);
size_t get_numberrows_bdd(PricerSolver *solver);

#ifdef __cplusplus
}
#endif

#endif // SOLVER_H
