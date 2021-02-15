#include <PricingStabilization.hpp>
#include <memory>
#include "PricerSolverArcTimeDP.hpp"
#include "PricerSolverBddBackward.hpp"
#include "PricerSolverBddForward.hpp"
#include "PricerSolverSimpleDP.hpp"
#include "PricerSolverZddBackward.hpp"
#include "PricerSolverZddForward.hpp"
#include "wctprivate.h"

extern "C" {
#include "scheduleset.h"

PricerSolverBase* copy_pricer_solver(PricerSolverBase* src,
                                     GPtrArray*        array,
                                     Parms*            parms) {
    switch (parms->pricing_solver) {
        case bdd_solver_simple:
            return new PricerSolverBddSimple(
                *dynamic_cast<PricerSolverBddSimple*>(src), array);
            break;
        case bdd_solver_cycle:
            return new PricerSolverBddCycle(
                *dynamic_cast<PricerSolverBddCycle*>(src), array);
            break;
        case zdd_solver_cycle:
            return new PricerSolverZddCycle(
                *dynamic_cast<PricerSolverZddCycle*>(src));
            break;
        case zdd_solver_simple:
            return new PricerSolverSimple(
                *dynamic_cast<PricerSolverSimple*>(src));
            break;
        case bdd_solver_backward_simple:
            return new PricerSolverBddBackwardSimple(
                *dynamic_cast<PricerSolverBddBackwardSimple*>(src), array);
            break;
        case bdd_solver_backward_cycle:
            return new PricerSolverBddBackwardCycle(
                *dynamic_cast<PricerSolverBddBackwardCycle*>(src), array);
            break;
        case zdd_solver_backward_simple:
            return new PricerSolverZddBackwardSimple(
                *dynamic_cast<PricerSolverZddBackwardSimple*>(src));
            break;
        case zdd_solver_backward_cycle:
            return new PricerSolverZddBackwardCycle(
                *dynamic_cast<PricerSolverZddBackwardCycle*>(src));
            break;
        case dp_solver:
            return new PricerSolverSimpleDp(
                *dynamic_cast<PricerSolverSimpleDp*>(src));
            break;
        case ati_solver:
            return new PricerSolverArcTimeDp(
                *dynamic_cast<PricerSolverArcTimeDp*>(src));
            break;
        case dp_bdd_solver:
            return new PricerSolverBddBackwardCycle(
                *dynamic_cast<PricerSolverBddBackwardCycle*>(src), array);
            break;
        default:
            return new PricerSolverBddCycle(
                *dynamic_cast<PricerSolverBddCycle*>(src), array);
    }
}

// void freeSolver(PricerSolver* src) {
//     delete src;
// }

// void deletePricerSolver(PricerSolver* solver) {
//     if (solver) {
//         delete solver;
//     }
// }

// int* get_take(PricerSolver* solver) {
//     return solver->get_take();
// }

// void iterate_zdd(PricerSolver* solver) {
//     solver->iterate_zdd();
// }

// int get_num_layers(PricerSolver* solver) {
//     return solver->get_num_layers();
// }

// size_t get_nb_vertices(PricerSolver* solver) {
//     return solver->get_nb_vertices();
// }

// size_t get_nb_edges(PricerSolver* solver) {
//     return solver->get_nb_edges();
// }

// void print_number_paths(PricerSolver* solver) {
//     solver->print_num_paths();
// }
// void print_dot_file(PricerSolver* solver, char* name) {
//     solver->create_dot_zdd(name);
// }

// void print_number_nodes_edges(PricerSolver* solver) {
//     solver->print_number_nodes_edges();
// }

// int reduce_cost_fixing(NodeData* pd) {
//     int    val = 0;
//     int    UB = pd->opt_sol->tw;
//     double LB = pd->LP_lower_bound_dual;

//     auto* aux_pi = static_cast<double*>(static_cast<void*>(pd->pi->data));
//     pd->solver->reduce_cost_fixing(aux_pi, UB, LB);
//     if (pd->depth == 0) {
//         pd->stat->size_graph_after_reduced_cost_fixing =
//             get_nb_edges(pd->solver);
//     }
//     return val;
// }

int NodeData::build_solve_mip() {
    int val = 0;

    solver->build_mip();

    return val;
}

int NodeData::construct_lp_sol_from_rmp() {
    int val = 0;

    val = lp_interface_get_nb_cols(RMP, &nb_cols);
    // CCcheck_val_2(val, "Failed to get nb cols");
    assert(nb_cols - id_pseudo_schedules == localColPool->len);

    lambda = CC_SAFE_REALLOC(lambda, nb_cols - id_pseudo_schedules, double);
    // CCcheck_NULL_2(pd->lambda, "Failed to allocate memory to pd->x");
    val = lp_interface_x(RMP, lambda, id_pseudo_schedules);
    // CCcheck_val_2(val, "Failed in lp_interface_x");
    solver->construct_lp_sol_from_rmp(lambda, localColPool, localColPool->len);

    return val;
}

int NodeData::generate_cuts() {
    // 1. add cuts to reformulation model
    int val = 0;

    val = solver->add_constraints();
    solver->insert_constraints_lp(this);
    solver->update_coeff_constraints();
    // 2. add cuts to lp relaxation wctlp
    // 3. adjust the pricing solver (add constraints to original model)
    return val;
}

int NodeData::delete_unused_rows_range(int first, int last) {
    int val = 0;

    lp_interface_deleterows(RMP, first, last);
    solver->remove_constraints(first, last - first + 1);
    solver_stab->remove_constraints(first, last - first + 1);
    g_array_remove_range(pi, first, last - first + 1);
    g_array_remove_range(rhs, first, last - first + 1);
    g_array_remove_range(lhs_coeff, first, last - first + 1);

    return val;
}

int NodeData::call_update_rows_coeff() {
    int val = 0;

    solver->update_rows_coeff(id_valid_cuts);

    return val;
}

int NodeData::check_schedule_set(ScheduleSet* set) {
    return static_cast<int>(solver->check_schedule_set(set->job_list));
}

void NodeData::make_schedule_set_feasible(ScheduleSet* set) {
    solver->make_schedule_set_feasible(set->job_list);
}

void NodeData::get_mip_statistics(enum MIP_Attr c) {
    Statistics* statistics = stat;
    switch (c) {
        case MIP_Attr_Nb_Vars:
            statistics->mip_nb_vars = solver->get_int_attr_model(c);
            break;
        case MIP_Attr_Nb_Constr:
            statistics->mip_nb_constr = solver->get_int_attr_model(c);
            break;
        case MIP_Attr_Obj_Bound:
            statistics->mip_obj_bound = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Obj_Bound_LP:
            statistics->mip_obj_bound_lp = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Mip_Gap:
            statistics->mip_rel_gap = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Run_Time:
            statistics->mip_run_time = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Status:
            statistics->mip_status = solver->get_int_attr_model(c);
            break;
        case MIP_Attr_Nb_Simplex_Iter:
            statistics->mip_nb_iter_simplex = solver->get_dbl_attr_model(c);
            break;
        case MIP_Attr_Nb_Nodes:
            statistics->mip_nb_nodes = solver->get_dbl_attr_model(c);
            break;
    }
}
}
