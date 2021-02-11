#include "parms.h"
#include <limits.h>
#include <string.h>
#include <util.h>
// #include "wctparms.h"

const double TIME_LIMIT = 7200.0;
const double ALPHA_STAB_INIT = 0.8;

parms::parms()
    : init_upper_bound(INT_MAX),
      bb_explore_strategy(min_bb_explore_strategy),
      use_strong_branching(min_strong_branching),
      bb_node_limit(0),
      nb_iterations_rvnd(3),
      branching_cpu_limit(TIME_LIMIT),
      alpha(ALPHA_STAB_INIT),
      pricing_solver(bdd_solver_backward_cycle),
      mip_solver(min_mip_solver),
      use_heuristic(min_use_heuristic),
      reduce_cost_fixing(min_reduced_cost),
      branchandbound(min_branch_and_bound),
      stab_technique(min_stab),
      print(min_print_size),
      jobfile(),
      pname(),
      nb_jobs(0),
      nb_machines(0) {}

int Parms::parms_set_file(std::string const& fname) {
    jobfile = std::string(fname);
    return 0;
}

int Parms::parms_set_pname(std::string const& fname) {
    pname = fname;
    return 0;
}

int Parms::parms_set_branching_cpu_limit(double limit) {
    branching_cpu_limit = limit;
    return 0;
}

int Parms::parms_set_branching_strategy(int strategy) {
    // bb_branch_strategy = strategy;
    return 0;
}

int Parms::parms_set_alpha(double alpha) {
    alpha = alpha;
    return 0;
}

int Parms::parms_set_mip_solver(int usage) {
    mip_solver = usage;
    return 0;
}

int Parms::parms_set_use_heuristic(int usage) {
    use_heuristic = usage;
    return 0;
}
int Parms::parms_set_reduce_cost(int usage) {
    reduce_cost_fixing = static_cast<reduced_cost_fixing_param>(usage);
    return 0;
}

int Parms::parms_set_strong_branching(int strong) {
    use_strong_branching = strong;
    return 0;
}

int Parms::parms_set_nb_machines(int _nb_machines) {
    nb_machines = _nb_machines;
    return 0;
}

int Parms::parms_set_nb_iterations_rvnd(int nb_iterations) {
    nb_iterations_rvnd = nb_iterations;
    return 0;
}

int Parms::parms_set_branchandbound(int bound) {
    branchandbound = bound;
    return 0;
}

int Parms::parms_set_bb_explore_strategy(int strategy) {
    bb_explore_strategy = static_cast<BBExploreStrategy>(strategy);
    return 0;
}

int Parms::parms_set_bb_node_limit(int node_limit) {
    bb_node_limit = node_limit;
    return 0;
}

int Parms::parms_set_stab_technique(int stab_technique) {
    stab_technique = static_cast<stab_techniques>(stab_technique);
    return 0;
}

int Parms::parms_set_print(int print) {
    print = print;
    return 0;
}

int Parms::parms_set_pricing_solver(int solver) {
    pricing_solver = solver;
    return 0;
}
