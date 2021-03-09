#include <fmt/core.h>
#include <memory>
#include "PricerSolverBase.hpp"
#include "scheduleset.h"
#include "wctprivate.h"

template <typename T = double>
int NodeData::construct_sol(OptimalSolution<T>* sol) {
    int                          val = 0;
    int                          nb_set = 1;
    std::shared_ptr<ScheduleSet> newset = std::make_shared<ScheduleSet>();

    newset->job_list = sol->jobs;
    // sol->jobs = nullptr;

    newset->total_weighted_completion_time = sol->cost;
    newset->total_processing_time = sol->C_max;
    newsets = newset;
    nb_new_sets = 1;

    return val;
}

int NodeData::solve_pricing() {
    int val = 0;

    update = 0;
    solver_stab->solve(LP_lower_bound_BB, pi.data(), lhs_coeff.data());

    if (solver_stab->get_update_stab_center()) {
        if (solver_stab->do_reduced_cost_fixing() &&
            parms->reduce_cost_fixing == yes_reduced_cost) {
            solver_stab->reduced_cost_fixing();
            check_schedules();
            delete_infeasible_schedules();
            solve_relaxation();
            double obj{};
            lp_interface_objval(RMP, &obj);
            solver_stab->update_continueLP(obj);
        }
    } else {
        if (!solver_stab->continueLP) {
            solver_stab->reduced_cost_fixing();
            check_schedules();
            delete_infeasible_schedules();
            solve_relaxation();
            double obj{};
            lp_interface_objval(RMP, &obj);
            solver_stab->update_continueLP(obj);
        }
    }

    if (solver_stab->get_reduced_cost() < -EPS_BOUND &&
        solver_stab->continueLP &&
        (solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND)) {
        auto sol = std::move(solver_stab->get_sol());
        val = construct_sol(&sol);
        CCcheck_val_2(val, "Failed in construction");
        val = add_lhs_scheduleset_to_rmp(newsets.get());
        newsets->id = localColPool.size();
        // g_ptr_array_add(localColPool, newsets);
        localColPool.emplace_back(newsets);
        // newsets = N;
        nb_new_sets = 0;
        nb_non_improvements = 0;
    } else {
        nb_new_sets = 0;
        nb_non_improvements++;
    }

CLEAN:
    return val;
}

int NodeData::solve_farkas_dbl() {
    int                     val = 0;
    OptimalSolution<double> s = solver->farkas_pricing(pi.data());
    update = 0;

    if (s.obj < EPS) {
        val = construct_sol(&s);
        lp_interface_write(RMP, "RMP.lp");

        // CCcheck_val_2(val, "Failed in constructing jobs");
        update = 1;
    } else {
        nb_new_sets = 0;
    }

    return val;
}
