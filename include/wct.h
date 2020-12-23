#ifndef _WCT_H
#define _WCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "defs.h"
#include "wctprivate.h"

/**
 * io.c
 */

int print_to_screen(Problem* problem);
int print_to_csv(Problem* problem);
int read_problem(Problem* problem);
int print_size_to_csv(Problem* problem, NodeData* pd);

/**
 * preprocess.c
 */

void calculate_Hmax(Problem* problem);
int  preprocess_data(Problem* problem);
int  find_division(Problem* problem);
void g_problem_summary_init(gpointer data, gpointer user_data);
void create_ordered_jobs_array(GPtrArray* a, GPtrArray* b);
void determine_jobs_order_interval(Problem* problem);

/**
 * greedy.c
 */

void update_bestschedule(Problem* problem, Solution* sol);

int construct_wspt(Job* jobarray, int nb_jobs, int nb_machines, Solution* sol);
int construct_feasible_solutions(Problem* problem);
int construct_edd(Problem* prob, Solution* sol);
int construct_spt(Problem* prob, Solution* sol);
int construct_random(Problem* prob, Solution* sol, GRand* rand_uniform);

int heuristic(Problem* prob);
int partlist_to_scheduleset(PartList*     part,
                            int           nb_part,
                            int           nb_jobs,
                            ScheduleSet** classes,
                            int*          column_count);

int prune_duplicated_sets(NodeData* pd);
/**
 * lowerbound.c
 */

int compute_lower_bound(NodeData* pd);
int compute_objective(NodeData* pd);
int print_x(NodeData* pd);
int calculate_x_e(NodeData* pd);
int calculate_nb_layers(NodeData* pd, int k);
int check_schedules(NodeData* pd);
int delete_infeasible_schedules(NodeData* pd);
int delete_old_schedules(NodeData* pd);
int delete_unused_rows(NodeData* pd);
int solve_relaxation(NodeData* pd);

void make_pi_feasible(NodeData* pd);
void make_pi_feasible_farkas_pricing(NodeData* pd);
int  add_newsets(NodeData* pd);

/** Help functions Glib */
void g_print_ages_col(gpointer data, gpointer user_data);
void g_grow_ages(gpointer data, gpointer user_data);
void g_make_pi_feasible(gpointer data, gpointer user_data);
void g_make_pi_feasible_farkas(gpointer data, gpointer user_data);

/**
 * model.c
 */

void g_add_col_to_lp(gpointer data, gpointer user_data);

int build_rmp(NodeData* pd);
int grab_integer_solution(NodeData* pd, double* x, double tolerance);
int add_scheduleset_to_rmp(ScheduleSet* set, NodeData* pd);
int add_lhs_scheduleset_to_rmp(ScheduleSet* set, NodeData* pd);
int get_solution_lp_lowerbound(NodeData* pd);
int add_artificial_var_to_rmp(NodeData* pd);

/**
 * wct.c
 */

int add_solution_to_colpool(Solution* sol, NodeData* pd);
int add_solution_to_colpool_and_lp(Solution* sol, NodeData* pd);

/**
 * solverwrapper.cc
 */

#ifdef __cplusplus
}
#endif

#endif
