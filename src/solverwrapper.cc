#include <wctprivate.h>
#include <PricerSolver.hpp>

extern "C" {

PricerSolverBase* newSolver(GPtrArray *jobs,
                            GPtrArray *ordered_jobs,
                            wctparms *parms) {
    switch (parms->pricing_solver) {
        case bdd_solver_simple:
            return new PricerSolverBddSimple(jobs, ordered_jobs);
        break;
        case bdd_solver_cycle:
            return new PricerSolverBddCycle(jobs, ordered_jobs);
        break;
        case zdd_solver_cycle:
            return new PricerSolverCycle(jobs, ordered_jobs);
        break;
        case zdd_solver_simple:
            return new PricerSolverZddSimple(jobs, ordered_jobs);
        break;
        case bdd_solver_backward_simple:
            return new PricerSolverBddBackwardSimple(jobs, ordered_jobs);
        break;
        case bdd_solver_backward_cycle:
            return new PricerSolverBddBackwardCycle (jobs, ordered_jobs);
        break;
        default:
            return new PricerSolverBddBackwardCycle(jobs, ordered_jobs);
    }
}

PricerSolverBase* newSolverDp(GPtrArray *_jobs, int _Hmax, wctparms *parms) {
    switch (parms->pricing_solver) {
        case dp_solver:
            return new PricerSolverSimpleDp(_jobs, _Hmax);
        break;
        case ati_solver:
            return new PricerSolverArcTimeDp(_jobs, _Hmax);
        default:
            return new PricerSolverSimpleDp(_jobs, _Hmax);
        break;
    }
}

void print_dot_file(PricerSolver *solver, char *name) {
    solver->create_dot_zdd(name);
}

void freeSolver(PricerSolver *src) { delete src; }

int evaluate_nodes(wctdata *pd) {
    int    val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound;
    int    nmachines = pd->problem->nmachines;

    pd->solver->evaluate_nodes(pd->pi, UB, LB, nmachines);

    return val;
}

int calculate_new_ordered_jobs(wctdata *pd) {
    int val = 0;

    pd->solver->calculate_new_ordered_jobs();

    return val;
}

int build_solve_mip(wctdata *pd) {
    int val = 0;

    // pd->solver->build_mip(pd->x_e);

    return val;
}

void print_number_nodes_edges(PricerSolver *solver) {
    solver->print_number_nodes_edges();
}

void deletePricerSolver(PricerSolver *solver) {
    if (solver) {
        delete solver;
    }
}

int calculate_table(PricerSolver *solver, wctparms *parms) {
    int val = 0;

    // solver->init_zdd_table();
    // solver->init_table_farkas();

    return val = 0;
}

void iterate_zdd(PricerSolver *solver) { solver->IterateZdd(); }

void print_number_paths(PricerSolver *solver) { solver->PrintNumberPaths(); }

size_t get_datasize(PricerSolver *solver) { return solver->get_datasize(); }

int get_nb_arcs_ati(PricerSolver *solver) { return solver->get_nb_arcs_ati(); }

size_t get_numberrows_zdd(PricerSolver *solver) { return solver->get_numberrows_zdd(); }

double get_edge_cost(PricerSolver *solver, int idx) { return solver->get_cost_edge(idx); }

int init_tables(PricerSolver *solver) {
    int val = 0;
    // solver->init_tables();
    return val;
}

void calculate_edges(PricerSolver *solver, scheduleset *set) {
    solver->calculate_edges(set);
}

void g_calculate_edges(gpointer data, gpointer user_data) {
    scheduleset *tmp = (scheduleset *) data;
    PricerSolver *solver = (PricerSolver *) user_data;

    solver->calculate_edges(tmp);
}

}