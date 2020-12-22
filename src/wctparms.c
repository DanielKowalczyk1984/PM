#include <limits.h>
#include <string.h>
#include <util.h>
#include <wctparms.h>

void parms_init(Parms* parms) {
    parms->init_upper_bound = INT_MAX;
    parms->bb_branch_strategy = min_bb_strategy;
    parms->bb_search_strategy = min_search_strategy;
    parms->bb_explore_strategy = min_bb_explore_strategy;
    parms->strong_branching = min_strong_branching;
    parms->mip_solver = min_mip_solver;
    parms->use_heuristic = min_use_heuristic;
    parms->reduce_cost_fixing = min_reduced_cost;
    parms->nb_iterations_rvnd = 3;
    parms->branchandbound = no;
    parms->stab_technique = no_stab;
    parms->print = min_print_size;
    parms->delete_edge_lists = 1;
    parms->jobfile = (char*)NULL;
    parms->pname = (char*)NULL;
    parms->upper_bounds_only = 0;
    parms->branching_cpu_limit = 3600.0;
    parms->alpha = 0.8;
    parms->pricing_solver = bdd_solver_backward_cycle;
}

void parms_free(Parms* parms) {
    CC_IFFREE(parms->jobfile, char);
    CC_IFFREE(parms->pname, char);
}

static int copy_string(char** dst, const char* src) {
    int val = 0;
    int len;
    len = (int)strlen(src) + 1;
    CC_IFFREE(*dst, char);
    *dst = (char*)CC_SAFE_MALLOC(len, char);
    CCcheck_NULL_2(*dst, "Failed to allocate dst");
    strcpy(*dst, src);
CLEAN:
    return val;
}

int parms_set_file(Parms* parms, const char* fname) {
    return copy_string(&(parms->jobfile), fname);
}

int parms_set_pname(Parms* parms, char* fname) {
    int val = copy_string(&(parms->pname), fname);
    CC_IFFREE(fname, char);
    return val;
}
/*Functions for setting the branching parameters*/
int parms_set_init_upper_bound(Parms* parms, int bound) {
    parms->init_upper_bound = bound;
    return 0;
}
int parms_set_branching_cpu_limit(Parms* parms, double limit) {
    parms->branching_cpu_limit = limit;
    return 0;
}

int parms_set_branching_strategy(Parms* parms, int strategy) {
    parms->bb_branch_strategy = strategy;
    return 0;
}

int parms_set_alpha(Parms* parms, double alpha) {
    parms->alpha = alpha;
    return 0;
}

int parms_set_search_strategy(Parms* parms, int strategy) {
    parms->bb_search_strategy = strategy;
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
    parms->mip_solver = usage;
    return 0;
}

int parms_set_strong_branching(Parms* parms, int strong) {
    parms->strong_branching = strong;
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

int parms_set_bb_strategy(Parms* parms, int strategy) {
    parms->bb_explore_strategy = strategy;
    return 0;
}

int parms_set_stab_technique(Parms* parms, int stab_technique) {
    parms->stab_technique = stab_technique;
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
