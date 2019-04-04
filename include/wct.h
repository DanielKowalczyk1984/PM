#ifndef _WCT_H
#define _WCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <wctprivate.h>
#include <defs.h>
#include <heap.h>

/**
 * io.c
 */

int print_to_screen(Problem *problem);
int print_to_csv(Problem *problem);
int read_problem(Problem *problem);
int print_size_to_csv(Problem *problem, wctdata *pd);

/**
 * preprocess.c
 */

void calculate_Hmax(Problem *problem);
int calculate_Hmin(
    int *durations, int nmachines, int njobs, int *perm, double *H);
int preprocess_data(Problem *problem);
int find_division(Problem *problem);
void g_problem_summary_init(gpointer data, gpointer user_data);
void create_ordered_jobs_array(GPtrArray *a, GPtrArray *b);
void determine_jobs_order_interval(Problem *problem);

/**
 * greedy.c
 */

void update_bestschedule(Problem *problem, Solution *sol);

int construct_wspt(Job *jobarray, int njobs, int nmachines, Solution *sol);
int construct_feasible_solutions(Problem *problem);
int construct_edd(Problem *prob, Solution *sol);
int construct_spt(Problem *prob, Solution *sol);
int construct_random(Problem *prob, Solution *sol, GRand *rand_uniform);

int heuristic_rpup(Problem *prob);
int partlist_to_scheduleset(
    partlist *part, int nbpart, int njobs, scheduleset **classes, int *ccount);

/**
 * branch_and_bound.c
 */

/*Help functions for branching*/
int insert_into_branching_heap(wctdata *pd, Problem *problem);
int skip_wctdata(wctdata *pd, Problem *problem);
int branching_msg(wctdata *pd, Problem *problem);
int branching_msg_cbfs(wctdata *pd, Problem *problem);
void free_elist(wctdata *cd, Parms *parms);
int prune_duplicated_sets(wctdata *pd);

/** Initialize BB tree */
void init_BB_tree(Problem *problem);

/** Conflict Branching functions */
int sequential_branching_conflict(Problem *problem);
/** Conflict Branching CBFS exploration */
int sequential_cbfs_branch_and_bound_conflict(Problem *problem);

/** help function for cbfs */
void insert_node_for_exploration(wctdata *pd, Problem *problem);
wctdata *get_next_node(Problem *problem);

int insert_frac_pairs_into_heap(wctdata *    pd,
                                       int *        nodepair_refs,
                                       double *     nodepair_weights,
                                       int          npairs,
                                       pmcheap *    heap);

/**
 * conflict_branching.c
 */
int create_branches_conflict(wctdata *pd, Problem *problem);

/**
 * lowerbound.c
 */

int compute_lower_bound(Problem *problem, wctdata *pd);
int compute_objective(wctdata *pd, Parms *parms);
int print_x(wctdata *pd);
int calculate_x_e(wctdata *pd);
int calculate_nblayers(wctdata *pd);

void make_pi_feasible(wctdata *pd);
void make_pi_feasible_farkas_pricing(wctdata *pd);

int add_newsets(wctdata *pd);

/** Help functions Glib */
void g_print_ages_col(gpointer data, gpointer user_data);
void g_grow_ages(gpointer data, gpointer user_data);
void g_make_pi_feasible(gpointer data, gpointer user_data);
void g_make_pi_feasible_farkas(gpointer data, gpointer user_data);

/**
 * model.c
 */

void g_add_col_to_lp(gpointer data, gpointer user_data);

int build_lp(wctdata *pd, int construct);
int grab_int_sol(wctdata *pd, double *x, double tolerance);
int addColToLP(scheduleset *set, wctdata *pd);
int get_solution_lp_lowerbound(wctdata *pd);


/**
 * wct.c
 */

void adapt_global_upper_bound(Problem *problem, int new_upper_bound);
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

int compute_schedule(Problem *problem);
int add_solution_to_colpool(Solution *sol, wctdata *pd);
int add_solution_to_colpool_and_lp(Solution *sol, wctdata *pd);


/**
 * solverwrapper.cc
 */

#ifdef __cplusplus
}
#endif

#endif
