#include <branch-and-bound/btree.h>
#include <fmt/core.h>
#include "PricerSolverBase.hpp"
#include "PricingStabilization.hpp"
#include "pricingstabilizationwrapper.h"
#include "scheduleset.h"
#include "stabilization.h"
#include "wct.h"
#include "wctparms.h"
#include "wctprivate.h"

template <typename T = double>
int construct_sol(NodeData* pd, OptimalSolution<T>* sol) {
    int          val = 0;
    int          nb_set = 1;
    ScheduleSet* newset = scheduleset_alloc_bis(pd->nb_jobs);
    CCcheck_NULL_3(newset, "Failed to allocate memory newset");

    newset->job_list = sol->jobs;
    sol->jobs = nullptr;

    newset->total_weighted_completion_time = sol->cost;
    newset->total_processing_time = sol->C_max;
    pd->newsets = newset;
    pd->nb_new_sets = 1;

CLEAN:
    if (val) {
        schedulesets_free(&(newset), &(nb_set));
    }

    return val;
}

extern "C" {

int solve_pricing(NodeData* pd) {
    int val = 0;

    pd->update = 0;
    Parms*  parms = pd->parms;
    double* pi = &g_array_index(pd->pi, double, 0);
    double* lhs = &g_array_index(pd->lhs_coeff, double, 0);

    pd->solver_stab->solve(pd->LP_lower_bound_BB, pi, lhs);

    if (call_get_update_stab_center(pd->solver_stab)) {
        if (call_do_reduced_fixing(pd->solver_stab) &&
            parms->reduce_cost_fixing == yes_reduced_cost) {
            pd->solver_stab->reduced_cost_fixing();
            check_schedules(pd);
            delete_infeasible_schedules(pd);
            solve_relaxation(pd);
            double obj{};
            lp_interface_objval(pd->RMP, &obj);
            call_update_continueLP(pd->solver_stab, obj);
        }
    } else {
        if (!call_get_continueLP(pd->solver_stab)) {
            pd->solver_stab->reduced_cost_fixing();
            check_schedules(pd);
            delete_infeasible_schedules(pd);
            solve_relaxation(pd);
            double obj{};
            lp_interface_objval(pd->RMP, &obj);
            call_update_continueLP(pd->solver_stab, obj);
        }
    }

    if (pd->solver_stab->get_reduced_cost() < -EPS_BOUND &&
        call_get_continueLP(pd->solver_stab) &&
        (call_get_eta_in(pd->solver_stab) <
         pd->upper_bound - 1.0 + EPS_BOUND)) {
        auto sol = std::move(pd->solver_stab->get_sol());
        val = construct_sol(pd, &sol);
        CCcheck_val_2(val, "Failed in construction");
        val = add_lhs_scheduleset_to_rmp(pd->newsets, pd);
        pd->newsets->id = pd->localColPool->len;
        g_ptr_array_add(pd->localColPool, pd->newsets);
        pd->newsets = NULL;
        pd->nb_new_sets = 0;
        pd->nb_non_improvements = 0;
    } else {
        pd->nb_new_sets = 0;
        pd->nb_non_improvements++;
    }

CLEAN:
    return val;
}

int solve_farkas_dbl(NodeData* pd) {
    int                     val = 0;
    OptimalSolution<double> s =
        pd->solver->farkas_pricing(&g_array_index(pd->pi, double, 0));
    pd->update = 0;

    if (s.obj < EPS) {
        val = construct_sol(pd, &s);
        lp_interface_write(pd->RMP, "RMP.lp");

        CCcheck_val_2(val, "Failed in constructing jobs");
        pd->update = 1;
    } else {
        pd->nb_new_sets = 0;
    }

CLEAN:
    return val;
}
}
