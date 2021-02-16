#include <fmt/core.h>
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>
#include "BranchBoundTree.hpp"
#include "Statistics.h"
#include "gurobi_c.h"
#include "scheduleset.h"
#include "solver.h"
#include "util.h"
#include "wctprivate.h"

void Problem::problem_init() {
    parms = Parms();
    stat = Statistics();
    statistics_init(&(stat));
    tree = nullptr;
    /** Job data */
    g_job_array = g_ptr_array_new_with_free_func(g_job_free);
    intervals = g_ptr_array_new_with_free_func(g_interval_free);
    opt_sol = nullptr;
    /** Job summary */
    nb_jobs = 0;
    p_sum = 0;
    pmax = 0;
    pmin = INT_MAX;
    dmax = INT_MIN;
    dmin = INT_MAX;
    off = 0;
    H_min = 0;
    H_max = INT_MAX;
    /*B&B info*/
    global_upper_bound = INT_MAX;
    global_lower_bound = 0;
    rel_error = DBL_MAX;
    root_lower_bound = 0;
    root_upper_bound = INT_MAX;
    root_rel_error = DBL_MAX;
    status = no_sol;
    // ColPool = g_ptr_array_new_with_free_func(g_scheduleset_free);
    root_pd = std::make_unique<NodeData>(this);
};

int debug = 0;

/*Information about debug*/
int dbg_lvl() {
    return debug;
}
void set_dbg_lvl(int dbglvl) {
    debug = dbglvl;
}

NodeData::NodeData(Problem* problem)
    : id(-1),
      depth(0),
      status(initialized),
      jobarray(nullptr),
      nb_jobs(problem->nb_jobs),
      nb_machines(problem->nb_machines),
      H_max(problem->H_max),
      H_min(problem->H_min),
      off(problem->off),
      ordered_jobs(g_ptr_array_new_with_free_func(g_free)),
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
      update(1),
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
      opt_sol(problem->opt_sol),
      pname("tmp") {
    /**init stab data */
}

NodeData::NodeData() {
    /*Initialization B&B data*/
    id = -1;
    depth = 0;
    status = initialized;
    pname = "temporary";
    /*Initialization node instance data*/
    nb_jobs = 0;
    H_max = 0;
    H_min = 0;
    off = 0;
    ordered_jobs = nullptr;
    jobarray = nullptr;
    /** Initialization data */
    upper_bound = INT_MAX;
    lower_bound = 0;
    LP_lower_bound = 0.0;
    LP_lower_min = DBL_MAX;
    rhs = std::vector<double>();
    /*Initialization  of the LP*/
    RMP = nullptr;
    lambda = std::vector<double>();
    // x_e = (double*)NULL;
    pi = std::vector<double>();
    slack = std::vector<double>();
    lhs_coeff = std::vector<double>();
    id_row = std::vector<int>();
    coeff_row = std::vector<double>();
    nb_rows = 0;
    nb_cols = 0;
    // init info cut generation
    max_nb_cuts = NB_CUTS;
    id_convex_constraint = 0;
    id_assignment_constraint = 0;
    id_valid_cuts = 0;
    id_art_var_assignment = 0;
    id_art_var_convex = 0;
    id_art_var_cuts = 0;
    id_pseudo_schedules = 0;

    /**init stab data */
    update = 1;
    iterations = 0;
    /*Initialization pricing_problem*/
    solver = nullptr;
    nb_non_improvements = 0;
    zero_count = 0;
    // bestcolors = (ScheduleSet*)NULL;
    best_schedule = std::vector<std::shared_ptr<ScheduleSet>>();
    // nb_best = 0;
    /**Column schedules */
    localColPool = std::vector<std::shared_ptr<ScheduleSet>>();
    column_status = std::vector<int>();
    /*Initialization max and retirement age*/
    maxiterations = NB_CG_ITERATIONS;
    retirementage = 0;
    /*initialization of branches*/
    branch_job = -1;
    parent = nullptr;
    completiontime = 0;
    less = -1;
    parms = nullptr;
}

std::unique_ptr<NodeData> NodeData::clone() {
    auto aux = std::make_unique<NodeData>();

    aux->parms = parms;
    aux->stat = stat;
    aux->opt_sol = opt_sol;
    aux->depth = depth + 1;

    /** Instance copy */
    aux->jobarray = jobarray;
    aux->nb_jobs = nb_jobs;
    aux->nb_machines = nb_machines;
    aux->H_max = H_max;
    aux->H_min = H_min;
    aux->off = off;

    /** copy info about intervals */
    aux->ordered_jobs =
        g_ptr_array_copy(ordered_jobs, g_copy_interval_pair, NULL);

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

    // aux->solver = copy_pricer_solver(solver, aux->ordered_jobs, aux->parms);
    aux->solver = solver->clone();
    // aux->solver->set_ordered_jobs(aux->ordered_jobs);
    // aux->ordered_jobs = aux->solver->get_ordered_jobs();

    aux->localColPool = localColPool;
    // g_ptr_array_copy(localColPool, g_copy_scheduleset, &(nb_jobs));

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

int set_id_and_name(NodeData* pd, int id, const char* fname) {
    int val = 0;
    pd->id = id;
    pd->pname = fmt::format("{}", fname);
    return val;
}

static void g_scheduleset_print_duplicated(gpointer data, gpointer user_data) {
    ScheduleSet* tmp = static_cast<ScheduleSet*>(data);
    GPtrArray*   tmp_array = tmp->job_list;

    fmt::print("TRANSSORT SET ");
}

void NodeData::prune_duplicated_sets() {
    auto equal_func = [](auto const& lhs, auto const& rhs) {
        return (*lhs) == (*rhs);
    };

    std::sort(localColPool.begin(), localColPool.end(),
              std::less<std::shared_ptr<ScheduleSet>>());
    localColPool.erase(
        std::unique(localColPool.begin(), localColPool.end(), equal_func),
        localColPool.end());

    int i = 0;
    for (auto& it : localColPool) {
        it->id = i++;
    }
}

void NodeData::add_solution_to_colpool(Solution* sol) {
    int val = 0;

    for (int i = 0; i < sol->nb_machines; ++i) {
        GPtrArray*                   machine = sol->part[i].machine;
        std::shared_ptr<ScheduleSet> tmp =
            std::make_shared<ScheduleSet>(machine);
        tmp->id = localColPool.size();
        // g_ptr_array_add(localColPool, tmp);
        localColPool.emplace_back(tmp);
    }
}

void NodeData::add_solution_to_colpool_and_lp(Solution* sol) {
    int val = 0;
    // ScheduleSet* tmp = nullptr;

    for (int i = 0; i < sol->nb_machines; ++i) {
        GPtrArray*                   machine = sol->part[i].machine;
        std::shared_ptr<ScheduleSet> tmp =
            std::make_shared<ScheduleSet>(machine);
        localColPool.emplace_back(tmp);
        // g_ptr_array_add(localColPool, tmp);
    }

    for (unsigned i = 0; i < localColPool.size(); ++i) {
        auto* tmp = localColPool[i].get();
        add_scheduleset_to_rmp(tmp);
    }
}
