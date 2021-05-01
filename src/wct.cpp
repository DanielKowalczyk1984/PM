#include <fmt/core.h>
#include <functional>
#include <limits>
#include <memory>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "lp.h"
#include "scheduleset.h"
#include "wctprivate.h"

NodeData::NodeData(Problem* problem)
    : depth(0UL),
      status(initialized),
      instance(problem->instance),
      nb_jobs(instance.nb_jobs),
      nb_machines(instance.nb_machines),
      RMP(lp_interface_create(nullptr), &lp_interface_delete),
      lambda(),
      pi(),
      slack(),
      rhs(),
      lhs_coeff(),
      id_row(),
      coeff_row(),
      nb_rows(0),
      nb_cols(0),
      max_nb_cuts(NB_CUTS),
      id_convex_constraint{},
      id_assignment_constraint{},
      id_valid_cuts{},
      id_art_var_convex{},
      id_art_var_assignment{},
      id_art_var_cuts{},
      id_next_var_cuts{},
      id_pseudo_schedules{},
      solver{},
      zero_count{},
      newsets(nullptr),
      nb_new_sets(0),
      column_status(),
      localColPool(),
      lower_bound(0),
      upper_bound(std::numeric_limits<int>::max()),
      LP_lower_bound(0.0),
      LP_lower_bound_dual(0.0),
      LP_lower_bound_BB(0.0),
      LP_lower_min(DBL_MAX),
      nb_non_improvements(0),
      iterations(0),
      solver_stab(nullptr),
      best_schedule(),
      best_objective(std::numeric_limits<int>::max()),
      retirementage(static_cast<int>(sqrt(nb_jobs)) + CLEANUP_ITERATION),
      branch_job(-1),
      completiontime(0),
      less(-1),
      parms(problem->parms),
      stat(problem->stat),
      opt_sol(problem->opt_sol),
      pname("tmp") {}

NodeData::NodeData(const NodeData& src)
    : depth(src.depth + 1),
      status(initialized),
      parms(src.parms),
      instance(src.instance),
      stat(src.stat),
      opt_sol(src.opt_sol),
      nb_jobs(src.nb_jobs),
      nb_machines(src.nb_machines),
      RMP(lp_interface_create(nullptr), &lp_interface_delete),
      lambda(),
      pi(),
      slack(),
      rhs(),
      lhs_coeff(),
      id_row(),
      coeff_row(),
      nb_rows(),
      nb_cols(),
      max_nb_cuts(src.max_nb_cuts),
      id_convex_constraint(src.id_convex_constraint),
      id_assignment_constraint(src.id_assignment_constraint),
      id_valid_cuts(src.id_valid_cuts),
      id_art_var_convex(src.id_art_var_convex),
      id_art_var_assignment(src.id_art_var_assignment),
      id_art_var_cuts(src.id_art_var_cuts),
      id_next_var_cuts(src.id_next_var_cuts),
      id_pseudo_schedules(src.id_pseudo_schedules),
      solver(src.solver->clone()),
      zero_count(),
      newsets(),
      nb_new_sets(),
      column_status(),
      localColPool(src.localColPool),
      lower_bound(src.lower_bound),
      upper_bound(src.upper_bound),
      LP_lower_bound(src.LP_lower_bound),
      LP_lower_bound_dual(src.LP_lower_bound_dual),
      LP_lower_bound_BB(src.LP_lower_bound_BB),
      LP_lower_min(src.LP_lower_min),
      nb_non_improvements(0),
      iterations(0),
      solver_stab(src.solver_stab->clone(solver.get(), pi)),
      best_schedule(src.best_schedule),
      best_objective(src.best_objective),
      retirementage(src.retirementage),
      branch_job(-1),
      completiontime(-1),
      less(-1) {}

std::unique_ptr<NodeData> NodeData::clone() const {
    auto aux = std::make_unique<NodeData>(*this);

    return aux;
}

void NodeData::prune_duplicated_sets() {
    auto equal_func = [](auto const& l, auto const& r) { return (*l) == (*r); };

    if (localColPool.empty() || localColPool.size() == 1) {
        return;
    }

    localColPool |=
        ranges::actions::sort(std::less<std::shared_ptr<ScheduleSet>>()) |
        ranges::actions::unique(equal_func);
    // std::ranges::sort(localColPool,
    // std::less<std::shared_ptr<ScheduleSet>>()); localColPool.erase(
    //     std::unique(localColPool.begin(), localColPool.end(), equal_func),
    //     localColPool.end());
}

void NodeData::add_solution_to_colpool(const Sol& sol) {
    for (auto& it : sol.machines) {
        localColPool.emplace_back(std::make_shared<ScheduleSet>(it));
    }
}
