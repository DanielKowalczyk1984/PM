#include <memory>
#include "PricerSolverBase.hpp"
#include "scheduleset.h"
#include "wctprivate.h"

template <typename T = double>
int NodeData::construct_sol(OptimalSolution<T>* sol) {
    int                          val = 0;
    std::shared_ptr<ScheduleSet> newset =
        std::make_shared<ScheduleSet>(std::move(*sol));

    newsets = newset;
    nb_new_sets = 1;

    return val;
}

int NodeData::solve_pricing() {
    int val = 0;

    solver_stab->solve(LP_lower_bound_BB, lhs_coeff.data());

    if (solver_stab->get_update_stab_center()) {
        if (solver_stab->do_reduced_cost_fixing() &&
            parms.reduce_cost_fixing == yes_reduced_cost) {
            solver_stab->reduced_cost_fixing();
            check_schedules();
            delete_infeasible_schedules();
            solve_relaxation();
            double obj{};
            lp_interface_objval(RMP.get(), &obj);
            solver_stab->update_continueLP(obj);
        }
    } else {
        if (!solver_stab->continueLP) {
            solver_stab->reduced_cost_fixing();
            check_schedules();
            delete_infeasible_schedules();
            solve_relaxation();
            double obj{};
            lp_interface_objval(RMP.get(), &obj);
            solver_stab->update_continueLP(obj);
        }
    }

    if (solver_stab->get_reduced_cost() < -EPS_BOUND &&
        solver_stab->continueLP &&
        (solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND)) {
        auto sol = std::move(solver_stab->get_sol());
        construct_sol(&sol);
        val = add_lhs_scheduleset_to_rmp(newsets.get());
        localColPool.emplace_back(newsets);
        nb_non_improvements = 0;
    } else {
        nb_non_improvements++;
    }

    return val;
}

void NodeData::solve_farkas_dbl() {
    OptimalSolution<double> s = solver->farkas_pricing(pi.data());

    if (s.obj < EPS) {
        construct_sol(&s);
        lp_interface_write(RMP.get(), "RMP.lp");
    } else {
        nb_new_sets = 0;
    }
}
