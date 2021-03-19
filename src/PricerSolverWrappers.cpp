#include "PricerSolverBase.hpp"
#include "scheduleset.h"
#include "wctprivate.h"

void NodeData::build_solve_mip() {
    solver->build_mip();
}

void NodeData::construct_lp_sol_from_rmp() {
    lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    assert(nb_cols - id_pseudo_schedules == localColPool.size());

    lambda.resize(nb_cols - id_pseudo_schedules, 0.0);
    lp_interface_x(RMP.get(), lambda.data(), id_pseudo_schedules);
    solver->construct_lp_sol_from_rmp(lambda.data(), localColPool,
                                      localColPool.size());
}

void NodeData::generate_cuts() {
    // 1. add cuts to reformulation model

    solver->add_constraints();
    solver->insert_constraints_lp(this);
    solver->update_coeff_constraints();
    // 2. add cuts to lp relaxation wctlp
    // 3. adjust the pricing solver (add constraints to original model)
}

int NodeData::delete_unused_rows_range(int first, int last) {
    int val = 0;

    lp_interface_deleterows(RMP.get(), first, last);
    solver->remove_constraints(first, last - first + 1);
    solver_stab->remove_constraints(first, last - first + 1);
    pi.erase(pi.begin(), pi.begin() + last - first + 1);
    pi.erase(rhs.begin(), rhs.begin() + last - first + 1);
    pi.erase(lhs_coeff.begin(), lhs_coeff.begin() + last - first + 1);

    return val;
}

int NodeData::call_update_rows_coeff() {
    int val = 0;

    solver->update_rows_coeff(id_valid_cuts);

    return val;
}

bool NodeData::check_schedule_set(ScheduleSet* set) {
    return solver->check_schedule_set(set->job_list);
}

void NodeData::make_schedule_set_feasible(ScheduleSet* set) {
    solver->make_schedule_set_feasible(set->job_list);
}
