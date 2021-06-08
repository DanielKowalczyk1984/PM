#include <fmt/core.h>
#include <memory>
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "scheduleset.h"
#include "wctprivate.h"

int NodeData::solve_pricing() {
    int val = 0;

    solver_stab->solve(LP_lower_bound_BB, lhs_coeff.data());

    if (solver_stab->get_update_stab_center()) {
        if (solver_stab->do_reduced_cost_fixing() &&
            parms.reduce_cost_fixing == yes_reduced_cost) {
            stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
            solver_stab->reduced_cost_fixing();
            stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
            // check_schedules();
            delete_infeasible_columns();
            solve_relaxation();
            double obj{};
            lp_interface_objval(RMP.get(), &obj);
            solver_stab->update_continueLP(obj);
        }
    } else {
        // if (!solver_stab->continueLP) {
        //     stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        //     solver_stab->reduced_cost_fixing();
        //     stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
        //     check_schedules();
        //     delete_infeasible_schedules();
        //     solve_relaxation();
        //     double obj{};
        //     lp_interface_objval(RMP.get(), &obj);
        //     solver_stab->update_continueLP(obj);
        // }
    }
    if (!solver_stab->continueLP) {
        stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        solver_stab->reduced_cost_fixing();
        stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
        delete_infeasible_columns();
        solve_relaxation();
        double obj{};
        lp_interface_objval(RMP.get(), &obj);
        solver_stab->update_continueLP(obj);
    }

    if (solver_stab->get_reduced_cost() < -EPS_BOUND &&
        solver_stab->continueLP &&
        (solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND)) {
        localColPool.emplace_back(
            std::make_shared<ScheduleSet>(std::move(solver_stab->get_sol())));
        val = add_lhs_scheduleset_to_rmp(localColPool.back().get());
        // nb_non_improvements = 0;
    } else {
        stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        solver_stab->reduced_cost_fixing();
        stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
        // check_schedules();
        delete_infeasible_columns();
        solve_relaxation();
        double obj{};
        lp_interface_objval(RMP.get(), &obj);
        solver_stab->update_continueLP(obj);
        // nb_non_improvements++;
    }

    return val;
}

void NodeData::solve_farkas_dbl() {
    OptimalSolution<double> s = solver->farkas_pricing(pi.data());

    if (s.obj < EPS) {
        localColPool.emplace_back(std::make_shared<ScheduleSet>(std::move(s)));
        add_lhs_scheduleset_to_rmp(localColPool.back().get());
    } else {
        nb_new_sets = 0;
    }
}
