#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H

#include "MIP_defs.hpp"
#include "Statistics.h"
#include "binomial-heap.h"
#include "interval.h"
#include "lp.h"
#include "pricingstabilizationwrapper.h"
#include "scheduleset.h"
#include "solver.h"
#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * wct data types nodes of branch and bound tree
 */
typedef enum {
    initialized = 0,
    LP_bound_estimated = 1,
    LP_bound_computed = 2,
    submitted_for_branching = 3,
    infeasible = 4,
    finished = 5,
} NodeDataStatus;

/**
 *  CONSTANTS NODEDATA STRUCTURE
 *
 */

#define NB_CUTS (2000)
#define NB_CG_ITERATIONS (1000000)
#define CLEANUP_ITERATION (30)
#define EPS (1e-6)
#define EPS_BOUND (1e-9)

/**
 * problem data
 */
typedef struct _Problem Problem;

/**
 * node data
 */
typedef struct _NodeData NodeData;

typedef struct BranchNodeBase  BranchNode;
typedef struct BranchBoundTree BranchBoundTree;

BranchBoundTree* new_branch_bound_tree(NodeData* data,
                                       int       _probtype,
                                       int       _isIntProb);

void delete_branch_bound_tree(BranchBoundTree* tree);
void call_branch_and_bound_explore(BranchBoundTree* tree);

struct _NodeData {
    // The id and depth of the node in the B&B tree
    int id;
    int depth;

    NodeDataStatus status;

    // The instance information
    GPtrArray* jobarray;
    int        nb_jobs;
    int        nb_machines;
    int        H_max;
    int        H_min;
    int        off;
    /** data about the intervals */
    GPtrArray* ordered_jobs;

    // The column generation lp information
    wctlp*  RMP;
    double* lambda;
    // double* x_e;

    GArray* pi;
    GArray* slack;
    GArray* rhs;
    GArray* lhs_coeff;
    GArray* id_row;
    GArray* coeff_row;

    int nb_rows;
    int nb_cols;

    // cut generation information
    int max_nb_cuts;
    int id_convex_constraint;
    int id_assignment_constraint;
    int id_valid_cuts;

    int id_art_var_convex;
    int id_art_var_assignment;
    int id_art_var_cuts;
    int id_next_var_cuts;
    int id_pseudo_schedules;

    // PricerSolver
    PricerSolver* solver;

    // Columns
    int          zero_count;
    ScheduleSet* newsets;
    int          nb_new_sets;
    int*         column_status;
    GPtrArray*   localColPool;

    int lower_bound;
    int upper_bound;

    double LP_lower_bound;
    double LP_lower_bound_dual;
    double LP_lower_bound_BB;
    double LP_lower_min;

    int nb_non_improvements;
    int iterations;

    /** Wentges smoothing technique */
    PricingStabilization* solver_stab;
    int                   update;

    // Best Solution
    GPtrArray* best_schedule;
    int        best_objective;

    // maxiterations and retireage
    int maxiterations;
    int retirementage;

    /** Branching strategies */
    int branch_job;
    int completiontime;
    int less;

    /**
     * ptr to the parent node
     */
    NodeData* parent;

    /** Some additional pointers to data needed */
    Parms*      parms;
    Statistics* stat;
    Solution*   opt_sol;

    char pname[MAX_PNAME_LEN];
};

/**
 * wct problem data type
 */

typedef enum {
    no_sol = 0,
    lp_feasible = 1,
    feasible = 2,
    meta_heuristic = 3,
    optimal = 4
} problem_status;

/*Initialization and free memory for the problem*/
// void problem_init(Problem* problem);
void problem_free(Problem* problem);

/*Initialize pmc data*/
void      nodedata_init(NodeData* pd, Problem* prob);
void      nodedata_init_null(NodeData* pd);
int       set_id_and_name(NodeData* pd, int id, const char* fname);
NodeData* new_node_data(NodeData* pd);

/*Free the Nodedata*/
void lp_node_data_free(NodeData* pd);
void nodedata_free(NodeData* pd);
void temporary_data_free(NodeData* pd);

/**
 * solver zdd
 */
int  build_solve_mip(NodeData* pd);
int  construct_lp_sol_from_rmp(NodeData* pd);
int  check_schedule_set(ScheduleSet* set, NodeData* pd);
void make_schedule_set_feasible(NodeData* pd, ScheduleSet* set);

void get_mip_statistics(NodeData* pd, enum MIP_Attr c);

/**
 * pricing algorithms
 */
int solve_pricing(NodeData* pd);
int solve_farkas_dbl(NodeData* pd);
int generate_cuts(NodeData* pd);
int delete_unused_rows_range(NodeData* pd, int first, int last);
int call_update_rows_coeff(NodeData* pd);
#ifdef __cplusplus
}
#endif
struct _Problem {
    Parms            parms;
    NodeData*        root_pd;
    BranchBoundTree* tree;
    /** Job data in EDD order */
    GPtrArray* g_job_array;
    GPtrArray* list_solutions;

    /** Summary of jobs */
    int nb_jobs;
    int p_sum;
    int pmax;
    int pmin;
    int dmax;
    int dmin;
    int H_min;
    int H_max;
    int off;
    int nb_machines;

    GPtrArray* intervals;

    int    nb_data_nodes;
    int    global_upper_bound;
    int    global_lower_bound;
    double rel_error;
    int    root_upper_bound;
    int    root_lower_bound;
    double root_rel_error;
    int    maxdepth;

    problem_status status;

    /* All partial schedules*/
    GPtrArray* ColPool;
    /* Best Solution*/
    Solution* opt_sol;
    /*Cpu time measurement + Statistics*/
    Statistics stat;

    void problem_init();
};

void _Problem::problem_init() {
    root_pd = CC_SAFE_MALLOC(1, NodeData);
    nodedata_init(root_pd, this);
    parms_init(&(parms));
    statistics_init(&(stat));
    tree = nullptr;
    /** Job data */
    g_job_array = g_ptr_array_new_with_free_func(g_job_free);
    intervals = g_ptr_array_new_with_free_func(g_interval_free);
    opt_sol = nullptr;
    /** Job summary */
    nb_jobs = 0;
    p_sum = 0;
    pmax = 0;
    pmin = INT_MAX;
    dmax = INT_MIN;
    dmin = INT_MAX;
    off = 0;
    H_min = 0;
    H_max = INT_MAX;
    /*B&B info*/
    nb_data_nodes = 0;
    global_upper_bound = INT_MAX;
    global_lower_bound = 0;
    rel_error = DBL_MAX;
    root_lower_bound = 0;
    root_upper_bound = INT_MAX;
    root_rel_error = DBL_MAX;
    status = no_sol;
    // br_heap_a = (BinomialHeap*)NULL;
    /*data of the problem*/
    set_id_and_name(root_pd, 0, "root_node");
    nb_data_nodes++;
    /*parms of the problem*/
    /** statistics of the problem */
    /*heap initialization*/
    /** initialize colPool */
    ColPool = g_ptr_array_new_with_free_func(g_scheduleset_free);
};

#endif
