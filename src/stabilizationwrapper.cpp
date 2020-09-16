#include "PricerSolverBase.hpp"
#include "PricingStabilization.hpp"
#include "fmt/core.h"
#include "scheduleset.h"
#include "stabilization.h"
#include "wctprivate.h"

template <typename T = double>
int construct_sol(NodeData* pd, OptimalSolution<T>* sol) {
    int          val = 0;
    int          nb_set = 1;
    ScheduleSet* newset = scheduleset_alloc_bis(pd->nb_jobs);
    CCcheck_NULL_3(newset, "Failed to allocate memory newset");

    for (unsigned i = 0; i < sol->jobs->len; ++i) {
        Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol->jobs, i));
        g_hash_table_add(newset->table, tmp_j);
    }
    newset->job_list = sol->jobs;
    sol->jobs = nullptr;
    newset->edge_list = nullptr;

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
    double* pi = &g_array_index(pd->pi, double, 0);
    double* lhs = &g_array_index(pd->lhs_coeff, double, 0);

    pd->solver_stab->solve(pd->eta_out, pi, lhs);

    if (pd->solver_stab->get_reduced_cost() < -1e-6) {
        pd->update = 1;
        auto sol = std::move(pd->solver_stab->get_sol());
        val = construct_sol(pd, &sol);
        CCcheck_val_2(val, "Failed in construction");
    } else {
        pd->nb_new_sets = 0;
    }

CLEAN:
    return val;
}

int solve_farkas_dbl(NodeData* pd) {
    int                     val = 0;
    OptimalSolution<double> s =
        pd->solver->farkas_pricing(&g_array_index(pd->pi, double, 0));
    pd->update = 0;

    if (s.obj < 1e-6) {
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
