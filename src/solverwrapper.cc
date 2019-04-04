#include <wctprivate.h>
#include <PricerSolver.hpp>

extern "C" {

PricerSolverBase* newSolver(GPtrArray *jobs,
                            int _num_machines,
                            GPtrArray *ordered_jobs,
                            wctparms *parms) {
    switch (parms->pricing_solver) {
        case bdd_solver_simple:
            return new PricerSolverBddSimple(jobs, _num_machines, ordered_jobs);
        break;
        case bdd_solver_cycle:
            return new PricerSolverBddCycle(jobs, _num_machines, ordered_jobs);
        break;
        // case zdd_solver_cycle:
        //     return new PricerSolverCycle(jobs, _num_machines, ordered_jobs);
        // break;
        // case zdd_solver_simple:
        //     return new PricerSolverZddSimple(jobs, _num_machines, ordered_jobs);
        // break;
        case bdd_solver_backward_simple:
            return new PricerSolverBddBackwardSimple(jobs, _num_machines, ordered_jobs);
        break;
        case bdd_solver_backward_cycle:
            return new PricerSolverBddBackwardCycle (jobs, _num_machines, ordered_jobs);
        break;
        default:
            return new PricerSolverBddBackwardCycle(jobs, _num_machines, ordered_jobs);
    }
}

PricerSolverBase* newSolverDp(GPtrArray *_jobs, int _num_machines, int _Hmax, wctparms *parms) {
    switch (parms->pricing_solver) {
        case dp_solver:
            return new PricerSolverSimpleDp(_jobs, _num_machines, _Hmax);
        break;
        case ati_solver:
            return new PricerSolverArcTimeDp(_jobs, _num_machines, _Hmax);
        default:
            return new PricerSolverSimpleDp(_jobs, _num_machines, _Hmax);
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

    pd->solver->evaluate_nodes(pd->pi, UB, LB);

    return val;
}

int calculate_new_ordered_jobs(wctdata *pd) {
    int val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound;

    pd->solver->calculate_new_ordered_jobs(pd->pi, UB, LB);

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

void iterate_zdd(PricerSolver *solver) { solver->iterate_zdd(); }

void print_number_paths(PricerSolver *solver) { solver->print_num_paths(); }

size_t get_datasize(PricerSolver *solver) { return solver->get_datasize(); }

int get_num_layers(PricerSolver *solver) { return solver->get_num_layers(); }

size_t get_size_graph(PricerSolver *solver) { return solver->get_size_graph(); }

double get_edge_cost(PricerSolver *solver, int idx) { return solver->get_cost_edge(idx); }

int init_tables(PricerSolver *solver) {
    int val = 0;
    // solver->init_tables();
    return val;
}

void calculate_edges(PricerSolver *solver, scheduleset *set) {
    solver->calculate_edges(set);
}

void construct_lp_sol_from_rmp(wctdata *pd) {
    pd->solver->construct_lp_sol_from_rmp(pd->x, pd->localColPool, pd->localColPool->len, pd->x_e);
}

void disjunctive_inequality(wctdata *pd, Solution *sol) {
    pd->solver->disjunctive_inequality(pd->x_e, sol);
}

void represent_solution(wctdata *pd, Solution *sol) {
    pd->solver->represent_solution(sol);
}

int check_schedule_set(scheduleset *set, wctdata *pd) {
    return pd->solver->check_schedule_set(set);
}

void g_calculate_edges(gpointer data, gpointer user_data) {
    scheduleset *tmp = (scheduleset *) data;
    PricerSolver *solver = (PricerSolver *) user_data;

    solver->calculate_edges(tmp);
}

}
