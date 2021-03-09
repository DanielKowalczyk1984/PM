#include <fmt/core.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <vector>
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "scheduleset.h"
#include "util.h"
#include "wctprivate.h"

NodeData::NodeData(Problem* problem)
    : id(-1),
      depth(0),
      status(initialized),
      instance((problem->instance)),
      nb_jobs(instance.nb_jobs),
      nb_machines(instance.nb_machines),
      RMP(nullptr),
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
      upper_bound(INT_MAX),
      LP_lower_bound(0.0),
      LP_lower_bound_dual(0.0),
      LP_lower_bound_BB(0.0),
      LP_lower_min(DBL_MAX),
      nb_non_improvements(0),
      iterations(0),
      solver_stab(nullptr),
      best_schedule(),
      best_objective(INT_MAX),
      maxiterations(NB_CG_ITERATIONS),
      retirementage(0),
      branch_job(-1),
      completiontime(0),
      less(-1),
      parent(nullptr),
      parms(&(problem->parms)),
      stat(&(problem->stat)),
      opt_sol(&(problem->opt_sol)),
      pname("tmp") {}

NodeData::NodeData(const Instance& _instance)
    : /*Initialization B&B data*/
      id{-1},
      depth{},
      status{initialized},
      instance(_instance),
      pname{"temporary"},
      /*Initialization node instance data*/
      nb_jobs{},
      nb_machines{},
      /** Initialization data */
      upper_bound{std::numeric_limits<int>::max()},
      lower_bound{},
      LP_lower_bound{},
      LP_lower_min{std::numeric_limits<double>::max()},
      rhs{},
      /*Initialization  of the LP*/
      RMP{},
      lambda{},
      pi{},
      slack{},
      lhs_coeff{},
      id_row{},
      coeff_row{},
      nb_rows{},
      nb_cols{},
      // init info cut generation
      max_nb_cuts{NB_CUTS},
      id_convex_constraint{},
      id_assignment_constraint{},
      id_valid_cuts{},
      id_art_var_assignment{},
      id_art_var_convex{},
      id_art_var_cuts{},
      id_pseudo_schedules{},

      /**init stab data */
      iterations{},
      /*Initialization pricing_problem*/
      solver{},
      nb_non_improvements{},
      zero_count{},
      best_schedule{},
      /**Column schedules */
      localColPool{},
      column_status{},
      /*Initialization max and retirement age*/
      maxiterations{NB_CG_ITERATIONS},
      retirementage{},
      /*initialization of branches*/
      branch_job{-1},
      parent{},
      completiontime{},
      less{-1},
      parms{} {}

std::unique_ptr<NodeData> NodeData::clone() const {
    auto aux = std::make_unique<NodeData>(instance);

    aux->parms = parms;
    aux->stat = stat;
    aux->opt_sol = opt_sol;
    aux->depth = depth + 1;

    /** Instance copy */
    aux->nb_jobs = nb_jobs;
    aux->nb_machines = nb_machines;

    /** Ids for model */
    aux->max_nb_cuts = max_nb_cuts;
    aux->id_convex_constraint = id_convex_constraint;
    aux->id_assignment_constraint = id_assignment_constraint;
    aux->id_valid_cuts = id_valid_cuts;

    aux->id_art_var_convex = id_art_var_convex;
    aux->id_art_var_assignment = id_art_var_assignment;
    aux->id_art_var_cuts = id_art_var_cuts;
    aux->id_next_var_cuts = id_next_var_cuts;
    aux->id_pseudo_schedules = id_pseudo_schedules;

    /** copy info about solver */

    aux->solver = solver->clone();

    aux->localColPool = localColPool;

    aux->lower_bound = lower_bound;
    aux->upper_bound = upper_bound;

    aux->LP_lower_bound = LP_lower_bound;
    aux->LP_lower_bound_dual = LP_lower_bound_dual;
    aux->LP_lower_bound_BB = LP_lower_bound_BB;
    aux->LP_lower_min = LP_lower_min;

    aux->iterations = 0;
    aux->nb_non_improvements = 0;

    aux->solver_stab = solver_stab->clone(aux->solver.get());
    /** copy info about best_schedule */
    aux->best_schedule = best_schedule;

    aux->maxiterations = maxiterations;
    aux->retirementage = retirementage;

    return aux;
}

void NodeData::lp_node_data_free() {
    /**
     * free all the gurobi data associated with the LP relaxation
     */
    if (RMP) {
        lp_interface_free(&(RMP));
    }

    nb_rows = 0;
    nb_cols = 0;
    max_nb_cuts = NB_CUTS;
    id_convex_constraint = 0;
    id_assignment_constraint = 0;
    id_valid_cuts = 0;
    id_art_var_assignment = 0;
    id_art_var_convex = 0;
    id_art_var_cuts = 0;
    id_pseudo_schedules = 0;
}

void NodeData::temporary_data_free() {
    lp_node_data_free();
}

void NodeData::prune_duplicated_sets() {
    auto equal_func = [](auto const& lhs, auto const& rhs) {
        return (*lhs) == (*rhs);
    };

    std::ranges::sort(localColPool, std::less<std::shared_ptr<ScheduleSet>>());
    localColPool.erase(
        std::unique(localColPool.begin(), localColPool.end(), equal_func),
        localColPool.end());

    int i = 0;
    for (auto& it : localColPool) {
        it->id = i++;
    }
}

void NodeData::add_solution_to_colpool(const Sol& sol) {
    for (auto& it : sol.machines) {
        std::shared_ptr<ScheduleSet> tmp = std::make_shared<ScheduleSet>(it);
        tmp->id = localColPool.size();
        localColPool.emplace_back(tmp);
    }
}

void NodeData::add_solution_to_colpool_and_lp(const Sol& sol) {
    add_solution_to_colpool(sol);

    for (auto& it : localColPool) {
        add_scheduleset_to_rmp(it.get());
    }
}
