#include "wctprivate.h"
#include <fmt/core.h>
#include <algorithm>
#include <boost/timer/timer.hpp>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>
#include "Instance.h"
#include "LocalSearch_new.h"
#include "PricerSolverArcTimeDP.hpp"
#include "PricerSolverBddBackward.hpp"
#include "PricerSolverBddForward.hpp"
#include "PricerSolverSimpleDP.hpp"
#include "PricerSolverZddBackward.hpp"
#include "PricerSolverZddForward.hpp"
#include "PricingStabilization.hpp"
#include "Solution.hpp"

Problem::~Problem() = default;

Problem::Problem(int argc, const char** argv)
    : parms(argc, argv),
      stat(),
      instance(&parms),
      tree(),
      root_pd(std::make_unique<NodeData>(this)),
      nb_jobs(instance.nb_jobs),
      nb_machines(instance.nb_machines),
      global_upper_bound(std::numeric_limits<int>::max()),
      global_lower_bound(),
      rel_error(std::numeric_limits<double>::max()),
      root_upper_bound(std::numeric_limits<int>::max()),
      root_lower_bound(),
      root_rel_error(std::numeric_limits<double>::max()),
      status(no_sol),
      opt_sol() {
    /**
     *@brief Finding heuristic solutions to the problem or start without
     *feasible solutions
     */
    if (parms.use_heuristic) {
        heuristic_new();
    } else {
        Sol best_sol(instance.nb_jobs, instance.nb_machines, instance.off);
        best_sol.construct_edd(instance.jobs);
        fmt::print("Solution Constructed with EDD heuristic:\n");
        best_sol.print_solution();
        best_sol.canonical_order(instance.intervals);
        fmt::print("Solution in canonical order: \n");
        best_sol.print_solution();
        opt_sol = best_sol;
    }
    /**
     * @brief Build DD at the root node
     *
     */
    CCutil_start_timer(&(stat.tot_build_dd));
    switch (parms.pricing_solver) {
        case bdd_solver_simple:
            root_pd->solver = std::make_unique<PricerSolverBddSimple>(instance);
            break;
        case bdd_solver_cycle:
            root_pd->solver = std::make_unique<PricerSolverBddCycle>(instance);
            break;
        case zdd_solver_cycle:
            root_pd->solver = std::make_unique<PricerSolverZddCycle>(instance);
            break;
        case zdd_solver_simple:
            root_pd->solver = std::make_unique<PricerSolverSimple>(instance);
            break;
        case bdd_solver_backward_simple:
            root_pd->solver =
                std::make_unique<PricerSolverBddBackwardSimple>(instance);
            break;
        case bdd_solver_backward_cycle:
            root_pd->solver =
                std::make_unique<PricerSolverBddBackwardCycle>(instance);
            break;
        case zdd_solver_backward_simple:
            root_pd->solver =
                std::make_unique<PricerSolverZddBackwardSimple>(instance);
            break;
        case zdd_solver_backward_cycle:
            root_pd->solver =
                std::make_unique<PricerSolverZddBackwardCycle>(instance);
        case dp_solver:
            root_pd->solver = std::make_unique<PricerSolverSimpleDp>(instance);
            break;
        case ati_solver:
            root_pd->solver = std::make_unique<PricerSolverArcTimeDp>(instance);
            break;
        case dp_bdd_solver:
            root_pd->solver = std::make_unique<PricerSolverSimpleDp>(instance);
            break;
        default:
            root_pd->solver =
                std::make_unique<PricerSolverBddBackwardCycle>(instance);
            break;
    }
    CCutil_stop_timer(&(stat.tot_build_dd), 0);
    stat.first_size_graph = root_pd->solver->get_nb_vertices();

    /**
     * @brief Initial stabilization method
     *
     */
    auto* tmp_solver = root_pd->solver.get();
    switch (parms.stab_technique) {
        case stab_wentgnes:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationStat>(tmp_solver);
            break;
        case stab_dynamic:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationDynamic>(tmp_solver);
            break;
        case stab_hybrid:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationHybrid>(tmp_solver);
            break;
        case no_stab:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationBase>(tmp_solver);
            break;
        default:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationStat>(tmp_solver);
            break;
    }
    root_pd->stat = &(stat);
    root_pd->opt_sol = &opt_sol;
    root_pd->solver->update_UB(opt_sol.tw);

    /**
     * @brief Initialization of the B&B tree
     *
     */
    tree = std::make_unique<BranchBoundTree>(std::move(root_pd), 0, 1);
    tree->explore();
}

NodeData::~NodeData() {
    temporary_data_free();
    // g_ptr_array_free(ordered_jobs, TRUE);
    // g_ptr_array_free(best_schedule, TRUE);
}

void Problem::solve() {
    tree->explore();
    if (parms.print) {
        print_to_csv();
    }
}

void Problem::heuristic_new() {
    auto                         ILS = nb_jobs / 2;
    auto                         IR = parms.nb_iterations_rvnd;
    boost::timer::auto_cpu_timer test;

    Sol best_sol(instance.nb_jobs, instance.nb_machines, instance.off);
    best_sol.construct_edd(instance.jobs);
    fmt::print("Solution Constructed with EDD heuristic:\n");
    best_sol.print_solution();
    best_sol.canonical_order(instance.intervals);
    fmt::print("Solution in canonical order: \n");
    best_sol.print_solution();

    /** Local Search */
    auto local = LocalSearchData(instance.nb_jobs, instance.nb_machines);
    local.RVND(best_sol);
    /** Perturbation operator */
    PerturbOperator perturb{};

    best_sol.canonical_order(instance.intervals);
    fmt::print("Solution after local search:\n");
    best_sol.print_solution();
    root_pd->add_solution_to_colpool(best_sol);

    Sol sol(instance.nb_jobs, instance.nb_machines, instance.off);
    int iterations = 0;
    for (auto i = 0UL; i < IR; ++i) {
        Sol sol1{instance.nb_jobs, instance.nb_machines, instance.off};
        sol1.construct_random_shuffle(instance.jobs);
        sol = sol1;

        for (auto j = 0UL; j < ILS; ++j) {
            local.RVND(sol1);
            ++iterations;
            if (sol1.tw < sol.tw) {
                sol = sol1;
                j = 0UL;
            }
            perturb(sol1);
        }

        if (sol.tw < best_sol.tw) {
            best_sol = sol;
            root_pd->add_solution_to_colpool(sol);
        }
    }

    fmt::print("Best new heuristics\n");
    best_sol.print_solution();
    opt_sol = best_sol;
}
