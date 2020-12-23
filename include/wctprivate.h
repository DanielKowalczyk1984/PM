#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H

#include "MIP_defs.hpp"
#include "Statistics.h"
#include "binomial-heap.h"
// #include "branch-and-boundwrapper.h"
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
 * problem data
 */
typedef struct _Problem Problem;

/**
 * node data
 */
typedef struct _NodeData NodeData;

typedef struct BranchNodeBase  BranchNode;
typedef struct BranchBoundTree BranchBoundTree;

BranchNode* new_branch_node(int _isRoot, NodeData* data);
void        delete_branch_node(BranchNode* node);
size_t      call_getDepth(BranchNode* state);
int         call_getDomClassID(BranchNode* state);
double      call_getObjValue(BranchNode* state);
double      call_getLB(BranchNode* state);
double      call_getUB(BranchNode* state);
int         call_getID(BranchNode* state);
int         call_getParentID(BranchNode* state);
void        call_setID(BranchNode* state, int i);
int         isDominated(BranchNode* state);
int         wasProcessed(BranchNode* state);

BranchBoundTree* new_branch_bound_tree(NodeData* data,
                                       int       _probtype,
                                       int       _isIntProb);
void             delete_branch_bound_tree(BranchBoundTree* tree);
void             call_branch_and_bound_explore(BranchBoundTree* tree);

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
    // int choose;
    // /** conflict */
    // int*      elist_same;
    // int       edge_count_same;
    // int*      elist_differ;
    // int       edge_count_differ;
    // NodeData* same_children;
    // int       nb_same;
    // NodeData* diff_children;
    // int       nb_diff;
    // Job *     v1, *v2;
    // /** ahv branching */
    // NodeData* duetime_child;
    // int       nb_duetime;
    // NodeData* releasetime_child;
    // int       nb_releasetime;
    // /** wide branching conflict */
    // int*       v1_wide;
    // int*       v2_wide;
    // int        nb_wide;
    // NodeData** same_children_wide;
    // NodeData** diff_children_wide;

    /** Branch info */
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
    /*heap variables*/
    BinomialHeap* br_heap_a;
    GPtrArray*    unexplored_states;
    GQueue*       non_empty_level_pqs;
    unsigned int  last_explored;
    int           mult_key;
    int           found;
    /*Cpu time measurement + Statistics*/
    Statistics stat;
};

/*Initialization and free memory for the problem*/
void problem_init(Problem* problem);
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

#endif
