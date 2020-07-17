#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include <scheduleset.h>
#include <wctparms.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PricerSolverBase PricerSolver;
PricerSolver*                   newSolver(GPtrArray* jobs, int _num_machines,
                                          GPtrArray* ordered_jobs, Parms* parms, int _Hmax,
                                          int* _take_jobs);
PricerSolver* newSolverDp(GPtrArray* _jobs, int _num_machines, int _Hmax,
                          Parms* parms);
PricerSolver* newSolverTIBdd(GPtrArray* _jobs, int _num_machines,
                             GPtrArray* _ordered_jobs, int* _take_jobs,
                             int _Hmax, Parms* parms);
void          copy_solver(PricerSolver** dest, PricerSolver* src);
void          freeSolver(PricerSolver* src);
void          deletePricerSolver(PricerSolver* solver);
int*          get_take(PricerSolver* solver);

void iterate_dd(PricerSolver* solver);
void iterate_zdd(PricerSolver* solver);

int    get_num_layers(PricerSolver* solver);
size_t get_nb_vertices(PricerSolver* solver);
size_t get_nb_edges(PricerSolver* solver);

void print_number_paths(PricerSolver* solver);
void print_dot_file(PricerSolver* solver, char* name);
void print_number_nodes_edges(PricerSolver* solver);

double get_edge_cost(PricerSolver* solver, int idx);

#ifdef __cplusplus
}
#endif
#endif  // INCLUDE_SOLVER_H_
