////////////////////////////////////////////////////////////////
//                                                            //
//  wctparms.h                                                //
//  wct                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef INCLUDE_WCTPARMS_H_
#define INCLUDE_WCTPARMS_H_

#ifdef __cplusplus
extern "C" {
#endif

const double TIME_LIMIT = 7200.0;
const double ALPHA_STAB_INIT = 0.8;

enum BBNodeSelection {
    min_search_strategy = 0,
    no_branching = min_search_strategy,
    min_lb_strategy = 1,
    dfs_strategy = 2,
    max_strategy = 3,
};

enum PricingSolver {
    bdd_solver_simple = 1,
    bdd_solver_cycle = 2,
    bdd_solver_backward_simple = 3,
    bdd_solver_backward_cycle = 4,
    zdd_solver_simple = 5,
    zdd_solver_cycle = 6,
    zdd_solver_backward_simple = 7,
    zdd_solver_backward_cycle = 8,
    dp_solver = 9,
    ati_solver = 10,
    dp_bdd_solver = 11,
};

enum BranchandBound {
    no = 0,
    yes = 1,
};

enum stab_techniques {
    no_stab = 0,
    stab_wentgnes = 1,
    stab_dynamic = 2,
    stab_hybrid = 3,
};

enum print {
    min_print_size = 0,
    use_print = 1,
};

enum BBBranchStrategy {
    min_bb_strategy = 0,
    conflict_strategy = min_bb_strategy,
    ahv_strategy = 1,
    cbfs_conflict_strategy = 2,
    cbfs_ahv_strategy = 3,
};

enum BBExploreStrategy {
    min_bb_explore_strategy = 0,
    bb_dfs_strategy = min_bb_explore_strategy,
    bb_bfs_strategy = 1,
    bb_brfs_strategy = 2,
    bb_cbfs_strategy = 3,
};

enum Strong_Branching {
    min_strong_branching = 0,
    use_strong_branching = min_strong_branching,
    no_strong_branching = 1,
};

enum MIP_solver {
    min_mip_solver = 0,
    no_mip_solver = min_mip_solver,
    use_mip_solver = 1,
};

enum reduced_cost_fixing_param {
    min_reduced_cost = 0,
    yes_reduced_cost = min_reduced_cost,
    no_reduced_cost = 1,
};

enum use_heuristic {
    min_use_heuristic = 0,
    yes_use_heuristic = min_use_heuristic,
    no_use_heuristic = 1,
};

typedef struct parms {
    /**
     * General parameters
     */
    int    init_upper_bound;
    int    bb_search_strategy;
    int    bb_explore_strategy;
    int    bb_branch_strategy;
    int    strong_branching;
    int    nb_iterations_rvnd;
    double branching_cpu_limit;
    double alpha;
    int    pricing_solver;
    int    mip_solver;
    int    use_heuristic;

    enum reduced_cost_fixing_param reduce_cost_fixing;

    /**
     * column generation
     */
    int branchandbound;
    int stab_technique;
    int print;

    int delete_edge_lists;

    char* jobfile;
    char* pname;

    int upper_bounds_only;
    int nb_jobs;
    int nb_machines;
} Parms;

/*Initialization and free memory*/
void parms_init(Parms* parms);
void parms_free(Parms* parms);

/*Functions for setting some parameters*/
int parms_set_init_upper_bound(Parms* parms, int bound);
int parms_set_branching_cpu_limit(Parms* parms, double limit);
int parms_set_alpha(Parms* parms, double alpha);
int parms_set_search_strategy(Parms* parms, int strategy);
int parms_set_branching_strategy(Parms* parms, int strategy);
int parms_set_strong_branching(Parms* parms, int strong);
int parms_set_mip_solver(Parms* parms, int usage);
int parms_set_use_heuristic(Parms* parms, int usage);
int parms_set_reduce_cost(Parms* parms, int usage);
int parms_set_nb_iterations_rvnd(Parms* parms, int nb_sol);

/**
 * Column Generation parameters
 */
int parms_set_pricing_solver(Parms* parms, int solver);
int parms_set_stab_technique(Parms* parms, int stab_technique);

/**
 * Branch-and-Bound parameters
 */
int parms_set_branchandbound(Parms* parms, int bound);
int parms_set_bb_strategy(Parms* parms, int strategy);
int parms_set_print(Parms* parms, int print);

/*Functions for defining the filesname*/
int parms_set_file(Parms* parms, const char* fname);
int parms_set_pname(Parms* parms, char* fname);
int parms_set_nb_machines(Parms* parms, int nb_machines);

#ifdef __cplusplus
}
#endif
#endif  // INCLUDE_WCTPARMS_H_
