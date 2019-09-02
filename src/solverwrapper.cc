#include <wctprivate.h>
#include "PricerSolverArcTimeDP.hpp"
#include "PricerSolverBddBackward.hpp"
#include "PricerSolverBddForward.hpp"
#include "PricerSolverSimpleDP.hpp"
#include "PricerSolverZddBackward.hpp"
#include "PricerSolverZddForward.hpp"

extern "C" {

PricerSolverBase* newSolver(GPtrArray* jobs, int _num_machines,
                            GPtrArray* ordered_jobs, wctparms* parms) {
    switch (parms->pricing_solver) {
        case bdd_solver_simple:
            return new PricerSolverBddSimple(jobs, _num_machines, ordered_jobs);
            break;
        case bdd_solver_cycle:
            return new PricerSolverBddCycle(jobs, _num_machines, ordered_jobs);
            break;
        case zdd_solver_cycle:
            return new PricerSolverZddCycle(jobs, _num_machines, ordered_jobs);
            break;
        case zdd_solver_simple:
            return new PricerSolverSimple(jobs, _num_machines, ordered_jobs);
            break;
        case bdd_solver_backward_simple:
            return new PricerSolverBddBackwardSimple(jobs, _num_machines,
                                                     ordered_jobs);
            break;
        case bdd_solver_backward_cycle:
            return new PricerSolverBddBackwardCycle(jobs, _num_machines,
                                                    ordered_jobs);
            break;
        case zdd_solver_backward_simple:
            return new PricerSolverZddBackwardSimple(jobs, _num_machines,
                                                     ordered_jobs);
            break;
        case zdd_solver_backward_cycle:
            return new PricerSolverZddBackwardCycle(jobs, _num_machines,
                                                    ordered_jobs);
        default:
            return new PricerSolverBddBackwardCycle(jobs, _num_machines,
                                                    ordered_jobs);
    }
}

PricerSolverBase* newSolverDp(GPtrArray* _jobs, int _num_machines, int _Hmax,
                              wctparms* parms) {
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

void print_dot_file(PricerSolver* solver, char* name) {
    solver->create_dot_zdd(name);
}

void freeSolver(PricerSolver* src) {
    delete src;
}

int evaluate_nodes(NodeData* pd) {
    int    val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound;

    pd->solver->evaluate_nodes(pd->pi, UB, LB);

    return val;
}

int reduce_cost_fixing(NodeData* pd) {
    int    val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound;

    pd->solver->reduce_cost_fixing(pd->pi, UB, LB);

    return val;
}

int build_solve_mip(NodeData* pd) {
    int val = 0;

    pd->solver->build_mip();

    return val;
}

void print_number_nodes_edges(PricerSolver* solver) {
    solver->print_number_nodes_edges();
}

void deletePricerSolver(PricerSolver* solver) {
    if (solver) {
        delete solver;
    }
}

void iterate_zdd(PricerSolver* solver) {
    solver->iterate_zdd();
}

void print_number_paths(PricerSolver* solver) {
    solver->print_num_paths();
}

int get_num_layers(PricerSolver* solver) {
    return solver->get_num_layers();
}

size_t get_size_graph(PricerSolver* solver) {
    return solver->get_size_graph();
}

void calculate_edges(PricerSolver* solver, ScheduleSet* set) {
    solver->calculate_edges(set);
}

void construct_lp_sol_from_rmp(NodeData* pd) {
    pd->solver->construct_lp_sol_from_rmp(pd->x, pd->localColPool,
                                          pd->localColPool->len, pd->x_e);
}

void disjunctive_inequality(NodeData* pd, Solution* sol) {
    pd->solver->disjunctive_inequality(pd->x_e, sol);
}

void represent_solution(NodeData* pd, Solution* sol) {
    pd->solver->represent_solution(sol);
}

int check_schedule_set(ScheduleSet* set, NodeData* pd) {
    return static_cast<int>(pd->solver->check_schedule_set(set->job_list));
}

void add_constraint(NodeData* pd, Job* job, int order) {
    pd->solver->add_constraint(job, pd->ordered_jobs, order);
}

void g_calculate_edges(gpointer data, gpointer user_data) {
    ScheduleSet*  tmp = (ScheduleSet*)data;
    PricerSolver* solver = (PricerSolver*)user_data;

    solver->calculate_edges(tmp);
}
}
