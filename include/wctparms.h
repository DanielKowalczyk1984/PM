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
    ati_solver = 10
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

enum BBSearchStrategy {
    min_bb_strategy = 0,
    conflict_strategy = min_bb_strategy,
    ahv_strategy = 1,
    cbfs_conflict_strategy = 2,
    cbfs_ahv_strategy = 3,
};

enum Strong_Branching {
    min_strong_branching = 0,
    use_strong_branching = min_strong_branching,
    no_strong_branching = 1,
};

typedef struct wctparms {
    /**
     * General parameters
     */
    int    init_upper_bound;
    int    bb_search_strategy;
    int    bb_branch_strategy;
    int    strong_branching;
    int    nb_iterations_rvnd;
    double branching_cpu_limit;
    double alpha;
    int pricing_solver;
    /**
     * scatter search
     */
    int    combine_method;
    int    scatter_search;
    double scatter_search_cpu_limit;
    /**
     * column generation
     */
    int branchandbound;
    int stab_technique;
    int print;

    int delete_elists;
    int delete_cclasses;

    char *jobfile;
    char *outfile;
    char *cclasses_infile;
    char *cclasses_outfile;
    char *color_infile;
    char *backupdir;

    int upper_bounds_only;
    int nmachines;
} Parms;

/*Initialization and free memory*/
void parms_init(Parms *parms);
void parms_free(Parms *parms);

/*Functions for setting some parameters*/
int parms_set_init_upper_bound(Parms *parms, int bound);
int parms_set_branching_cpu_limit(Parms *parms, double limit);
int parms_set_alpha(Parms *parms, double alpha);
int parms_set_search_strategy(Parms *parms, int strategy);
int parms_set_branching_strategy(Parms *parms, int strategy);
int parms_set_strong_branching(Parms *parms, int strong);
int parms_set_nb_iterations_rvnd(Parms *parms, int nb_sol);

/**
 * column generation
 */
int parms_set_branchandbound(Parms *parms, int bound);
int parms_set_stab_technique(Parms *parms, int stab_technique);
int parms_set_print(Parms *parms, int print);
int parms_set_pricing_solver(Parms *parms, int solver);

/**
 * scatter search
 */
int parms_set_scatter_search(Parms *parms, int scatter);
int parms_set_combine_method(Parms *parms, int combine_method);
int parms_set_scatter_search_cpu_limit(Parms *parms, double limit);

/*Functions for defining the filesname*/
int parms_set_outfile(Parms *parms, const char *fname);
int parms_set_file(Parms *parms, const char *fname);
int parms_set_backupdir(Parms *parms, const char *fname);
int parms_set_nmachines(Parms *parms, int nmachines);

#ifdef __cplusplus
}
#endif
#endif  // INCLUDE_WCTPARMS_H_
