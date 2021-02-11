#ifndef _WCT_H
#define _WCT_H

#include "defs.h"
#include "wctprivate.h"
int parse_cmd(int argc, const char** argv, Parms* parms);

#ifdef __cplusplus
extern "C" {
#endif

/**
@brief Parse the commands with docopt (see parse_cmd.cpp for definition)
 *
 * @param argc
 * @param argv
 * @param parms
 * @return int
 */

/**
 * io.c
 */

int print_to_screen(Problem* problem);
int print_to_csv(Problem* problem);
int read_problem(Problem* problem);

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
