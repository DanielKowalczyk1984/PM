#include "parms.h"
#include <limits.h>
#include <string.h>
#include <util.h>
// #include "wctparms.h"

const double TIME_LIMIT = 7200.0;
const double ALPHA_STAB_INIT = 0.8;

// void parms_init(Parms* parms) {
//     parms->init_upper_bound = INT_MAX;
//     parms->bb_branch_strategy = min_bb_strategy;
//     parms->bb_explore_strategy = min_bb_explore_strategy;
//     parms->bb_node_limit = 0;
//     parms->use_strong_branching = min_strong_branching;
//     parms->mip_solver = min_mip_solver;
//     parms->use_heuristic = min_use_heuristic;
//     parms->reduce_cost_fixing = min_reduced_cost;
//     parms->nb_iterations_rvnd = 3;
//     parms->branchandbound = min_branch_and_bound;
//     parms->stab_technique = min_stab;
//     parms->print = min_print_size;
//     parms->jobfile = (char*)NULL;
//     parms->pname = (char*)NULL;
//     parms->branching_cpu_limit = TIME_LIMIT;
//     parms->alpha = ALPHA_STAB_INIT;
//     parms->pricing_solver = bdd_solver_backward_cycle;
// }

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
      jobfile(nullptr),
      pname(nullptr),
      nb_jobs(0),
      nb_machines(0) {}

parms::~parms() {
    CC_IFFREE(jobfile, char);
    CC_IFFREE(pname, char);
}

// void parms_free(Parms* parms) {}

static int copy_string(char** dst, const char* src) {
    int val = 0;
    int len = 0;
    len = (int)strlen(src) + 1;
    CC_IFFREE(*dst, char);
    *dst = CC_SAFE_MALLOC(len, char);
    CCcheck_NULL_2(*dst, "Failed to allocate dst");
    strcpy(*dst, src);
CLEAN:
    return val;
}

int parms_set_file(Parms* parms, const char* fname) {
    return copy_string(&(parms->jobfile), fname);
}

int parms_set_pname(Parms* parms, const char* fname) {
    int val = copy_string(&(parms->pname), fname);
    return val;
}

int parms_set_branching_cpu_limit(Parms* parms, double limit) {
    parms->branching_cpu_limit = limit;
    return 0;
}

int parms_set_branching_strategy(Parms* parms, int strategy) {
    // parms->bb_branch_strategy = strategy;
    return 0;
}

int parms_set_alpha(Parms* parms, double alpha) {
    parms->alpha = alpha;
    return 0;
}

int parms_set_mip_solver(Parms* parms, int usage) {
    parms->mip_solver = usage;
    return 0;
}

int parms_set_use_heuristic(Parms* parms, int usage) {
    parms->use_heuristic = usage;
    return 0;
}
int parms_set_reduce_cost(Parms* parms, int usage) {
    parms->reduce_cost_fixing = static_cast<reduced_cost_fixing_param>(usage);
    return 0;
}

int parms_set_strong_branching(Parms* parms, int strong) {
    parms->use_strong_branching = strong;
    return 0;
}

int parms_set_nb_machines(Parms* parms, int nb_machines) {
    parms->nb_machines = nb_machines;
    return 0;
}

int parms_set_nb_iterations_rvnd(Parms* parms, int nb_iterations) {
    parms->nb_iterations_rvnd = nb_iterations;
    return 0;
}

int parms_set_branchandbound(Parms* parms, int bound) {
    parms->branchandbound = bound;
    return 0;
}

int parms_set_bb_explore_strategy(Parms* parms, int strategy) {
    parms->bb_explore_strategy = static_cast<BBExploreStrategy>(strategy);
    return 0;
}

int parms_set_bb_node_limit(Parms* parms, int node_limit) {
    parms->bb_node_limit = node_limit;
    return 0;
}

int parms_set_stab_technique(Parms* parms, int stab_technique) {
    parms->stab_technique = static_cast<stab_techniques>(stab_technique);
    return 0;
}

int parms_set_print(Parms* parms, int print) {
    parms->print = print;
    return 0;
}

int parms_set_pricing_solver(Parms* parms, int solver) {
    parms->pricing_solver = solver;
    return 0;
}
