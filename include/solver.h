#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include <wctparms.h>
#include <solution.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PricerSolverBase PricerSolver;
PricerSolver *newSolver(GPtrArray *jobs, GPtrArray *ordered_jobs, wctparms *parms);
PricerSolver *newSolverDp(GPtrArray *_jobs, int _Hmax, wctparms *parms);
PricerSolver *newSolverDP(GPtrArray *interval_list, int njobs, int **sum_p);
// PricerSolver *copySolver(PricerSolver *src);
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
void print_number_paths(PricerSolver *solver);

void set_release_due_time(PricerSolver *solver, Job *jobarray);

size_t get_datasize(PricerSolver *solver);
int get_nb_arcs_ati(PricerSolver *solver);
size_t get_numberrows_zdd(PricerSolver *solver);
size_t get_numberrows_bdd(PricerSolver *solver);

#ifdef __cplusplus
}
#endif
#endif  // INCLUDE_SOLVER_H_
