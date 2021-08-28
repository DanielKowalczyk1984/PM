#include "NodeData.h"  // for NodeData
#include <fmt/core.h>
#include <array>                                // for array
#include <cmath>                                // for sqrt
#include <concepts/concepts.hpp>                // for return_t
#include <cstddef>                              // for size_t
#include <functional>                           // for less, function
#include <limits>                               // for numeric_limits
#include <memory>                               // for shared_ptr, unique_ptr
#include <range/v3/action/action.hpp>           // for action_closure, opera...
#include <range/v3/action/sort.hpp>             // for sort, sort_fn
#include <range/v3/action/unique.hpp>           // for unique, unique_fn
#include <range/v3/functional/bind_back.hpp>    // for bind_back_fn_
#include <range/v3/functional/comparisons.hpp>  // for equal_to
#include <range/v3/functional/compose.hpp>      // for composed
#include <utility>                              // for move
#include <vector>                               // for vector
#include "Column.h"                             // for ScheduleSet
#include "Instance.h"                           // for Instance
#include "Parms.h"                              // for Parms
#include "PricerSolverBase.hpp"                 // for PricerSolverBase
#include "PricingStabilization.hpp"             // for PricingStabilizationBase
#include "Problem.h"                            // for Problem
#include "Solution.hpp"                         // for Sol
#include "lp.h"                                 // for lp_interface_create

NodeData::NodeData(Problem* problem)
    : depth(0UL),
      status(initialized),
      parms(problem->parms),
      instance(problem->instance),
      stat(problem->stat),
      opt_sol(problem->opt_sol),
      pname("tmp"),
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
      LP_lower_min(std::numeric_limits<double>::max()),
      nb_non_improvements(0),
      iterations(0UL),
      solver_stab(nullptr),
      retirementage(static_cast<int>(sqrt(static_cast<double>(nb_jobs))) +
                    CLEANUP_ITERATION),
      branch_job(),
      completiontime(0),
      less(-1) {}

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
      iterations(0UL),
      solver_stab(src.solver_stab->clone(solver.get(), pi)),
      retirementage(src.retirementage),
      branch_job(),
      completiontime(-1),
      less(-1) {}

std::unique_ptr<NodeData> NodeData::clone() const {
    auto aux = std::make_unique<NodeData>(*this);

    return aux;
}

NodeData::~NodeData() = default;

std::unique_ptr<NodeData> NodeData::clone(size_t _j, int _t, bool _left) const {
    auto aux = std::make_unique<NodeData>(*this);
    aux->branch_job = _j;
    aux->completiontime = _t;
    aux->less = _left;
    aux->solver->split_job_time(_j, _t, _left);
    aux->build_rmp();
    aux->delete_infeasible_columns();
    aux->solve_relaxation();
    aux->estimate_lower_bound(2 * instance.nb_machines);
    return aux;
}

std::array<std::unique_ptr<NodeData>, 2> NodeData::create_child_nodes(size_t _j,
                                                                      int _t) {
    return std::array<std::unique_ptr<NodeData>, 2>{clone(_j, _t, false),
                                                    clone(_j, _t, true)};
}

std::array<std::unique_ptr<NodeData>, 2> NodeData::create_child_nodes(size_t _j,
                                                                      long _t) {
    auto aux_t = static_cast<int>(_t);
    return std::array<std::unique_ptr<NodeData>, 2>{clone(_j, aux_t, false),
                                                    clone(_j, aux_t, true)};
}

void NodeData::prune_duplicated_sets() {
    // auto equal_func = [](auto const& l, auto const& r) { return (*l) == (*r);
    // };

    if (localColPool.empty() || localColPool.size() == 1) {
        return;
    }

    localColPool |=
        ranges::actions::sort(std::less<>{},
                              [](auto& tmp) { return tmp->job_list; }) |
        ranges::actions::unique(ranges::equal_to{},
                                [](auto& tmp) { return tmp->job_list; });
    // std::ranges::sort(localColPool,
    // std::less<std::shared_ptr<ScheduleSet>>()); localColPool.erase(
    //     std::unique(localColPool.begin(), localColPool.end(), equal_func),
    //     localColPool.end());
}

void NodeData::add_solution_to_colpool(const Sol& sol) {
    for (auto& it : sol.machines) {
        localColPool.emplace_back(std::make_shared<Column>(it));
    }
}

double NodeData::get_score_value() {
    switch (parms.scoring_value) {
        case (size_scoring_value):
            return static_cast<double>(solver->get_nb_edges());
            break;
        case (nb_paths_scoring_value):
            return static_cast<double>(solver->print_num_paths());
            break;
        default:
            return LP_lower_bound;
    }
}
