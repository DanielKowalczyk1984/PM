#include "PricerSolverArcTimeDP.hpp"
#include "PricerSolverBddBackward.hpp"
#include "PricerSolverBddForward.hpp"
#include "PricerSolverSimpleDP.hpp"
#include "PricerSolverZddBackward.hpp"
#include "PricerSolverZddForward.hpp"
#include "wctprivate.h"

extern "C" {
#include "scheduleset.h"

PricerSolverBase* newSolver(GPtrArray* jobs,
                            int        _num_machines,
                            GPtrArray* ordered_jobs,
                            parms*     parms,
                            int        _Hmax,
                            int*       _take_jobs,
                            double     _UB) {
    switch (parms->pricing_solver) {
        case bdd_solver_simple:
            return new PricerSolverBddSimple(jobs, _num_machines, ordered_jobs,
                                             parms->pname, _Hmax, _take_jobs,
                                             _UB);
            break;
        case bdd_solver_cycle:
            return new PricerSolverBddCycle(jobs, _num_machines, ordered_jobs,
                                            parms->pname, _Hmax, _take_jobs,
                                            _UB);
            break;
        case zdd_solver_cycle:
            return new PricerSolverZddCycle(jobs, _num_machines, ordered_jobs,
                                            parms->pname, _UB);
            break;
        case zdd_solver_simple:
            return new PricerSolverSimple(jobs, _num_machines, ordered_jobs,
                                          parms->pname, _UB);
            break;
        case bdd_solver_backward_simple:
            return new PricerSolverBddBackwardSimple(jobs, _num_machines,
                                                     ordered_jobs, parms->pname,
                                                     _Hmax, _take_jobs, _UB);
            break;
        case bdd_solver_backward_cycle:
            return new PricerSolverBddBackwardCycle(jobs, _num_machines,
                                                    ordered_jobs, parms->pname,
                                                    _Hmax, _take_jobs, _UB);
            break;
        case zdd_solver_backward_simple:
            return new PricerSolverZddBackwardSimple(
                jobs, _num_machines, ordered_jobs, parms->pname, _UB);
            break;
        case zdd_solver_backward_cycle:
            return new PricerSolverZddBackwardCycle(
                jobs, _num_machines, ordered_jobs, parms->pname, _UB);
        default:
            return new PricerSolverBddBackwardCycle(jobs, _num_machines,
                                                    ordered_jobs, parms->pname,
                                                    _Hmax, _take_jobs, _UB);
    }
}

PricerSolverBase* newSolverDp(GPtrArray* _jobs,
                              int        _num_machines,
                              int        _Hmax,
                              parms*     parms,
                              double     _UB) {
    switch (parms->pricing_solver) {
        case dp_solver:
            return new PricerSolverSimpleDp(_jobs, _num_machines, _Hmax,
                                            parms->pname, _UB);
            break;
        case ati_solver:
            return new PricerSolverArcTimeDp(_jobs, _num_machines, _Hmax,
                                             parms->pname, _UB);
        case dp_bdd_solver:
            return new PricerSolverSimpleDp(_jobs, _num_machines, _Hmax,
                                            parms->pname, _UB);
        default:
            return new PricerSolverSimpleDp(_jobs, _num_machines, _Hmax,
                                            parms->pname, _UB);
            break;
    }
}

void freeSolver(PricerSolver* src) {
    delete src;
}

void deletePricerSolver(PricerSolver* solver) {
    if (solver) {
        delete solver;
    }
}

int* get_take(PricerSolver* solver) {
    return solver->get_take();
}

void iterate_zdd(PricerSolver* solver) {
    solver->iterate_zdd();
}

int get_num_layers(PricerSolver* solver) {
    return solver->get_num_layers();
}

size_t get_nb_vertices(PricerSolver* solver) {
    return solver->get_nb_vertices();
}

size_t get_nb_edges(PricerSolver* solver) {
    return solver->get_nb_edges();
}

void print_number_paths(PricerSolver* solver) {
    solver->print_num_paths();
}
void print_dot_file(PricerSolver* solver, char* name) {
    solver->create_dot_zdd(name);
}

void print_number_nodes_edges(PricerSolver* solver) {
    solver->print_number_nodes_edges();
}

int evaluate_nodes(NodeData* pd) {
    int    val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound;

    pd->solver->evaluate_nodes(&g_array_index(pd->pi, double, 0), UB, LB);

    return val;
}

int reduce_cost_fixing(NodeData* pd) {
    int    val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound_dual;

    pd->solver->reduce_cost_fixing(&g_array_index(pd->pi, double, 0), UB, LB);
    pd->problem->size_graph_after_reduced_cost_fixing =
        get_nb_edges(pd->solver);
    return val;
}

int build_solve_mip(NodeData* pd) {
    int val = 0;

    pd->solver->build_mip();

    return val;
}

int construct_lp_sol_from_rmp(NodeData* pd) {
    int val = 0;
    int nb_cols;

    val = lp_interface_get_nb_cols(pd->RMP, &nb_cols);
    CCcheck_val_2(val, "Failed to get nb cols");
    assert(nb_cols - pd->id_pseudo_schedules == pd->localColPool->len);

    pd->lambda =
        CC_SAFE_REALLOC(pd->lambda, nb_cols - pd->id_pseudo_schedules, double);
    CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
    val = lp_interface_x(pd->RMP, pd->lambda, pd->id_pseudo_schedules);
    CCcheck_val_2(val, "Failed in lp_interface_x");
    pd->solver->construct_lp_sol_from_rmp(pd->lambda, pd->localColPool,
                                          pd->localColPool->len);

CLEAN:
    return val;
}

void generate_cuts(NodeData* pd) {
    // 1. add cuts to reformulation model
    PricerSolver* pricing_solver = pd->solver;
    pricing_solver->add_constraints();
    pricing_solver->insert_constraints_lp(pd);
    pricing_solver->update_coeff_constraints();
    // 2. add cuts to lp relaxation wctlp
    // 3. adjust the pricing solver (add constraints to original model)
}

void represent_solution(NodeData* pd, Solution* sol) {
    pd->solver->represent_solution(sol);
}

int check_schedule_set(ScheduleSet* set, NodeData* pd) {
    return static_cast<int>(pd->solver->check_schedule_set(set->job_list));
}

void make_schedule_set_feasible(NodeData* pd, ScheduleSet* set) {
    PricerSolver* solver = pd->solver;
    solver->make_schedule_set_feasible(set->job_list);
}

void add_constraint(NodeData* pd, Job* job, int order) {
    pd->solver->add_constraint(job, pd->ordered_jobs, order);
}

void get_mip_statistics(NodeData* pd, enum MIP_Attr c) {
    Problem*          problem = pd->problem;
    PricerSolverBase* solver = pd->solver;
    switch (c) {
        case MIP_Attr_Nb_Vars:
            problem->mip_nb_vars = solver->get_int_attr_model(c);
            break;
        case MIP_Attr_Nb_Constr:
            problem->mip_nb_constr = solver->get_int_attr_model(c);
            break;
        case MIP_Attr_Obj_Bound:
            problem->mip_obj_bound = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Obj_Bound_LP:
            problem->mip_obj_bound_lp = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Mip_Gap:
            problem->mip_rel_gap = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Run_Time:
            problem->mip_run_time = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Status:
            problem->mip_status = solver->get_int_attr_model(c);
            break;
        case MIP_Attr_Nb_Simplex_Iter:
            problem->mip_nb_iter_simplex = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Nb_Nodes:
            problem->mip_nb_nodes = solver->get_dbl_attr_model(c);
            break;
    }
}
}
