#ifndef _WCT_H
#define _WCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <defs.h>
#include <wctprivate.h>
#include <heap.h>

/**
 * io.c
 */

int print_to_screen(wctproblem *problem);
int print_to_csv(wctproblem *problem);
int read_problem(wctproblem *problem);
int print_size_to_csv(wctproblem *problem, wctdata *pd);

/**
 * preprocess.c
 */

int calculate_Hmax(Job *jobarray, int nmachines, int njobs);
int calculate_Hmin(
    int *durations, int nmachines, int njobs, int *perm, double *H);
int preprocess_data(wctproblem *problem);
int find_division(wctproblem *problem);
void g_problem_summary_init(gpointer data, gpointer user_data);

/**
 * greedy.c
 */

void update_bestschedule(wctproblem *problem, solution *sol);

int construct_wspt(Job *jobarray, int njobs, int nmachines, solution *sol);
int construct_feasible_solutions(wctproblem *problem);
int construct_edd(wctproblem *prob, solution *sol);
int construct_spt(wctproblem *prob, solution *sol);
int construct_random(wctproblem *prob, solution *sol, GRand *rand_uniform);

int heuristic_rpup(wctproblem *prob);
int partlist_to_scheduleset(
    partlist *part, int nbpart, int njobs, scheduleset **classes, int *ccount);

/**
 * branch_and_bound.c
 */

/*Help functions for branching*/
int insert_into_branching_heap(wctdata *pd, wctproblem *problem);
int skip_wctdata(wctdata *pd, wctproblem *problem);
int branching_msg(wctdata *pd, wctproblem *problem);
int branching_msg_cbfs(wctdata *pd, wctproblem *problem);
void free_elist(wctdata *cd, wctparms *parms);
int prune_duplicated_sets(wctdata *pd);

/** Initialize BB tree */
void init_BB_tree(wctproblem *problem);

/** Conflict Branching functions */
int sequential_branching_conflict(wctproblem *problem);
/** Conflict Branching CBFS exploration */
int sequential_cbfs_branch_and_bound_conflict(wctproblem *problem);

/** help function for cbfs */
void insert_node_for_exploration(wctdata *pd, wctproblem *problem);
wctdata *get_next_node(wctproblem *problem);

int insert_frac_pairs_into_heap(wctdata *    pd,
                                       const double x[],
                                       int *        nodepair_refs,
                                       double *     nodepair_weights,
                                       int          npairs,
                                       pmcheap *    heap);

/**
 * conflict_branching.c
 */
int create_branches_conflict(wctdata *pd, wctproblem *problem);

/**
 * lowerbound.c
 */

int compute_lower_bound(wctproblem *problem, wctdata *pd);
int compute_objective(wctdata *pd, wctparms *parms);

void make_pi_feasible(wctdata *pd);
void make_pi_feasible_farkas_pricing(wctdata *pd);

int add_newsets(wctdata *pd);

/**
 * model.c
 */

int build_lp(wctdata *pd, int construct);
int grab_int_sol(wctdata *pd, double *x, double tolerance);


/**
 * wct.c
 */

void adapt_global_upper_bound(wctproblem *problem, int new_upper_bound);
/** compute row-index v1 and column-index v2 from array-index.*/
static inline void inodepair_ref_key(int *v1, int *v2, int index) {
    *v2 = (int)floor(sqrt(2 * ((double)index) + 0.25) - 0.5);
    *v1 = index - (*v2 * (*v2 + 1) / 2);
}

/** help functions for heap srong branching */
static inline int nodepair_ref_key(int v1, int v2) {
    /* We store only the elements of the upper right triangle within the
     vcount x vcount matrix. */
    assert(v1 <= v2);
    return v2 * (v2 + 1) / 2 + v1;
}

int compute_schedule(wctproblem *problem);


/**
 * wide_branching.c
 */

int create_branches_wide(wctdata *pd, wctproblem *problem);
int sequential_branching_wide(wctproblem *problem);
int branching_msg_wide(wctdata *pd, wctproblem *problem);

/**
 * solverwrapper.cc
 */

/**
 * Stabilization techniques
 */
int solve_stab(wctdata *pd, wctparms *parms);
int solve_stab_dynamic(wctdata *pd, wctparms *parms);

/**
 * solver zdd
 */
int solve_dynamic_programming_ahv(wctdata *pd);
int solve_weight_dbl_bdd(wctdata *pd);
int solve_weight_dbl_zdd(wctdata *pd);
int solve_pricing(wctdata *pd, wctparms *parms);
int solve_farkas_dbl(wctdata *pd);
int solve_farkas_dbl_DP(wctdata *pd);
void print_dot_file(PricerSolver *solver, char *name);

#ifdef __cplusplus
}
#endif

#endif
