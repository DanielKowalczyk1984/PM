#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include <wctparms.h>
#include <solution.h>
#include <scheduleset.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PricerSolverBase PricerSolver;
PricerSolver *newSolver(GPtrArray *jobs, int _num_machines, GPtrArray *ordered_jobs, wctparms *parms);
PricerSolver *newSolverDp(GPtrArray *_jobs, int _num_machines, int _Hmax, wctparms *parms);
void copy_solver(PricerSolver **dest, PricerSolver *src);
void freeSolver(PricerSolver *src);
void deletePricerSolver(PricerSolver *solver);

// PricerSolver *copySolver(PricerSolver *src);
int init_tables(PricerSolver *solver);
int calculate_table(PricerSolver *solver, wctparms *parms);



void iterate_dd(PricerSolver *solver);
void iterate_zdd(PricerSolver *solver);

// int get_nb_arcs_ati(PricerSolver *solver);
// size_t get_datasize(PricerSolver *solver);
int get_num_layers(PricerSolver *solver);
size_t get_size_graph(PricerSolver *solver);
// size_t get_numberrows_zdd(PricerSolver *solver);
// size_t get_numberrows_bdd(PricerSolver *solver);

void print_number_paths(PricerSolver *solver);
void print_dot_file(PricerSolver *solver, char *name);
void print_number_nodes_edges(PricerSolver *solver);

void calculate_edges(PricerSolver *solver, scheduleset *set);
void g_calculate_edges(gpointer data, gpointer user_data);
double get_edge_cost(PricerSolver *solver, int idx);

#ifdef __cplusplus
}
#endif
#endif  // INCLUDE_SOLVER_H_
