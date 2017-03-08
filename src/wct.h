#ifndef _WCT_H
#define _WCT_H


#include "wctprivate.h"
#include "defs.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * greedy.c
 */

struct _pair_job_machine {
    int job;
    int machine;
} ;

typedef struct _pair_job_machine pair_job_machine;


int random_rcl_assignment(Job *jobarray, int njobs, int nmachines,
                          solution *new_sol, GRand *rand_);
int construct_wspt(Job *jobarray, int njobs, int  nmachines,
                   solution  *new_sol);
int random_assignment(Job *jobarray, int njobs, int nmachines,
                      solution *new_sol, GRand *rand_);

void update_bestschedule(wctproblem *problem, solution *new_sol);
int construct_feasible_solutions(wctproblem *problem);
int construct_edd(wctproblem *prob, solution *sol);
int construct_spt(wctproblem *prob, solution *sol);
int heuristic_rpup(wctproblem *prob);

/*
compare_functions
 */
gint compare_func1(gconstpointer a, gconstpointer b, void *user_data);


/**
 * lowerbound.c
 */

typedef struct _MACHINE {
    double totweight;
    double totcompletion;
} MACHINE;


int lowerbound_cw(Job *array, int njobs, int nmachines);
int lowerbound_cp(Job *array, int njobs, int nmachines);
int lowerbound_eei(Job *array, int njobs, int nmachines);
int lowerbound_ak(Job *array, int njobs, int nmachines);

/**
 * compare_functions
 */

int compare_cw(BinomialHeapValue a, BinomialHeapValue b);


/**
 * Compare functions
 */

int order_refset(const void *a, const void *b);
int order_makespan(const void *a, const void *b, void *data);
int order_distance(const void *a, const void *b);
int order_makespan_list(const void *a, const void *b);


/**
 * For each functions
 */

void distance_min_max(void *data, void *user_data);
void max_dist(void *data, void *user_data);
void free_sol(void *data, void *user_data);
void assign_iter(void *data, void *user_data);
void print_sol(void *data, void *user_data);
void print_makespan(void *data, void *user_data);
void refset_dist(void *data, void *user_data);
void for_each_comp_fitness(void *data, void *user_data);

/**
 * Print functions
 */

void print_pool(SS *scatter_search);
void print_refset(SS *scatter_search);
void print_distance(SS *scatter_search);
void print_pool_n(SS *scatter_search, int n);
void print_pool_makespan(SS *scatter_search);
void print_refset_makespan(SS *scatter_search);
void print_list1(SS *scatter_search);

/**
 * localsearch.c
 */

/* local search methods */
typedef struct _slope_t
{
    int b1; /* begin of segment*/
    int b2; /* end of segment*/
    int c;  /* value of the function at b1*/
    int alpha; /* slope of the function */
} slope_t;

typedef struct _processing_list_data{
  int pos;
  int p;
} processing_list_data;

typedef struct _local_search_data
{
  int nmachines;
  int njobs;
  int **W;
  GList *** g;
  processing_list_data ** processing_list_f;
  processing_list_data ** processing_list_b;
  processing_list_data ** processing_list_inter1;
  processing_list_data ** processing_list_inter2;
  int updated;
} local_search_data;

local_search_data *local_search_data_init(solution *sol);
void local_search_data_free(local_search_data *data);

/** Preperation of the data */
int local_search_create_processing_list(solution *sol, local_search_data *data,
                                        int l) ;
int local_search_create_processing_list_b(solution *sol, local_search_data *data,
                                        int l) ;
int local_search_create_W(solution *sol, local_search_data *data);
int local_search_create_g(solution *sol, local_search_data *data);
/** Search for the best intra insertion */
void local_search_forward_insertion(solution *sol, local_search_data *data, int l);
void local_search_backward_insertion(solution *sol, local_search_data *data,int l) ;
int local_search_create_processing_list_insertion_inter(solution *sol, local_search_data *data,
                                        int l) ;
/** update of the intra insertion */
void local_seach_update_insertion(solution *sol, int i_best, int j_best, int k_best, int l, int improvement);
void local_search_swap_intra(solution *sol, local_search_data *data,
                                    int l1, int l2);
void local_search_insertion_inter(solution *sol, local_search_data *data,
                                    int l) ;
void local_search_update_insertion_inter(solution *sol, int i_best,
        int j_best, int k_best, int kk_best, int l, int improvement);
void local_search_swap_inter(solution *sol, local_search_data *data,
                                    int l1, int l2);


/**
 * wct.c
 */

/*Initialization and free memory for the problem*/
void wctproblem_init(wctproblem *problem);
void wctproblem_free(wctproblem *problem);

int compute_schedule(wctproblem *problem);

/*Computation functions*/
int compute_lower_bound(wctproblem *problem, wctdata *pd);
int sequential_branching(wctproblem *problem);
int create_branches(wctdata *pd, wctproblem *problem);
int check_integrality(wctdata *pd, int nmachine, int *result);
int build_lp(wctdata *pd, int construct);
void make_pi_feasible(wctdata *pd);
void make_pi_feasible_farkas_pricing(wctdata *pd);
int heur_colors_with_stable_sets(wctdata *pd);
int compute_objective(wctdata *pd, wctparms *parms);
/** Wide Branching functions */
int create_branches_wide(wctdata *pd, wctproblem *problem);
int sequential_branching_wide(wctproblem *problem);
/** Conflict Branching functions */
int create_branches_conflict(wctdata *pd, wctproblem *problem);
int sequential_branching_conflict(wctproblem *problem);
/** Conflict Branching functions */
int create_branches_ahv(wctdata *pd, wctproblem *problem);
int sequential_branching_ahv(wctproblem *problem);
/** Conflict Branching CBFS exploration */
int sequential_cbfs_branch_and_bound_conflict(wctproblem *problem);
/** AHV branching CBFS exploration */
int sequential_cbfs_branch_and_bound_ahv(wctproblem *problem);
/** Initialize BB tree */
void init_BB_tree(wctproblem *problem);


/*Help functions for branching*/
int insert_into_branching_heap(wctdata *pd, wctproblem *problem);
int skip_wctdata(wctdata *pd, wctproblem *problem);
int branching_msg(wctdata *pd, wctproblem *problem);
int branching_msg_wide(wctdata *pd, wctproblem *problem);
int branching_msg_cbfs(wctdata *pd, wctproblem *problem);
int branching_msg_cbfs_ahv(wctdata *pd, wctproblem *problem);
int collect_same_children(wctdata *pd);
int collect_diff_children(wctdata *pd);
void temporary_data_free(wctdata *pd);
int new_eindex(int v, int v1, int v2);
int prune_duplicated_sets(wctdata *pd);
int double2int(int *kpc_pi, int *scalef, const double *pi, int vcount);
double safe_lower_dbl(int numerator, int denominator);
int add_newsets(wctdata *pd);
int add_to_init(wctproblem *problem, Scheduleset *sets, int nbsets);
int delete_to_bigcclasses(wctdata *pd, int capacity);
int add_some_maximal_stable_sets(wctdata *pd, int number);
int insert_into_branching_heap(wctdata *pd, wctproblem *problem);
/** help function for cbfs */
void insert_node_for_exploration(wctdata *pd, wctproblem *problem);
wctdata *get_next_node(wctproblem *problem);

int remove_finished_subtreebis(wctdata *child);

void make_pi_feasible(wctdata *pd);

/*Backup*/
int backup_wctdata(wctdata *pd, wctproblem *problem);

/**
 * wct.c
 */

int read_problem(wctproblem *problem);

/** Preprocess data */
gint compare_readytime(gconstpointer a, gconstpointer b);
int calculate_ready_due_times(Job *jobarray, int njobs, int nmachines,
                              double Hmin);
int calculate_Hmax(Job *jobarray, int nmachines, int njobs);
int calculate_Hmin(int *durations, int nmachines, int njobs, int *perm,
                   double *H);
int preprocess_data(wctproblem *problem);


/**
 * heurdiving.c
 */

int heur_exec(wctproblem *problem, wctdata *pd, int *result);
void heur_init(wctdata *pd);
void heur_free(wctdata *pd);

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
int solvedblzdd(wctdata *pd);
int solvedblbdd(wctdata *pd);
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
