#include <ext/alloc_traits.h>        // for __alloc_traits<>::value_type
#include <memory>                    // for shared_ptr, unique_ptr, make_shared
#include <utility>                   // for move
#include <vector>                    // for vector
#include "Column.h"                  // for ScheduleSet
#include "NodeData.h"                // for NodeData
#include "OptimalSolution.hpp"       // for OptimalSolution
#include "Parms.h"                   // for Parms, yes_reduced_cost
#include "PricerSolverBase.hpp"      // for PricerSolverBase
#include "PricingStabilization.hpp"  // for PricingStabilizationBase
#include "Statistics.h"              // for Statistics, Statistics::reduced_...
#include "lp.h"                      // for lp_interface_objval

int NodeData::solve_pricing() {
    int val = 0;

    solver_stab->solve(LP_lower_bound_BB, lhs_coeff.data());

    if (solver_stab->get_update_stab_center()) {
        if (solver_stab->do_reduced_cost_fixing() && parms.reduce_cost_fixing) {
            stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
            if (solver_stab->reduced_cost_fixing()) {
                delete_infeasible_columns();
            }
            stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
            // solve_relaxation();
            // double obj{};
            // lp_interface_objval(RMP.get(), &obj);
            // solver_stab->update_continueLP(obj);
        }
    }

    if (!solver_stab->continueLP) {
        stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        if (solver_stab->reduced_cost_fixing()) {
            delete_infeasible_columns();
        }
        stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
        // solve_relaxation();
        // double obj{};
        // lp_interface_objval(RMP.get(), &obj);
        // solver_stab->update_continueLP(obj);
    }

    if (solver_stab->get_reduced_cost() < -EPS_BOUND &&
        solver_stab->continueLP &&
        (solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND)) {
        localColPool.emplace_back(
            std::make_shared<Column>(std::move(solver_stab->get_sol())));
        val = add_lhs_scheduleset_to_rmp(localColPool.back().get());
        // solve_relaxation();
        // double obj{};
        // lp_interface_objval(RMP.get(), &obj);
        // solver_stab->update_continueLP(obj);
        // nb_non_improvements = 0;
    } else {
        stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        if (solver_stab->reduced_cost_fixing()) {
            delete_infeasible_columns();
        }
        stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
        // solve_relaxation();
        // double obj{};
        // lp_interface_objval(RMP.get(), &obj);
        // solver_stab->update_continueLP(obj);
    }
    solve_relaxation();
    // double obj{};
    // lp_interface_objval(RMP.get(), &obj);
    solver_stab->update_continueLP(LP_lower_bound);

    return val;
}

void NodeData::solve_farkas_dbl() {
    OptimalSolution<double> s = solver->farkas_pricing(pi.data());

    if (s.obj < EPS) {
        localColPool.emplace_back(std::make_shared<Column>(std::move(s)));
        add_lhs_scheduleset_to_rmp(localColPool.back().get());
    } else {
        nb_new_sets = 0;
    }
}
