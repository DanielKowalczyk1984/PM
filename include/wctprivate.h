#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H

#include <Scheduleset.h>
#include <binomial-heap.h>
#include <lp.h>
#include <solution.h>
#include <util.h>
#include <wctparms.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Solver data types
 */

typedef struct PricerSolver PricerSolver;
PricerSolver *newSolver(Job *jobarray, int nbjobs, int Hmin, int Hmax);
PricerSolver *newSolverDP(Job *jobarray, int nbjobs, int Hmin, int Hmax);
PricerSolver *copySolver(PricerSolver *src);
void freeSolver(PricerSolver *src);
int calculate_table(PricerSolver *solver, wctparms *parms);
void deletePricerSolver(PricerSolver *solver);
int add_conflict_constraints(PricerSolver *solver,
                             wctparms *    parms,
                             int *         elist_same,
                             int           ecount_same,
                             int *         elist_differ,
                             int           ecount_differ);
int free_conflict_constraints(PricerSolver *solver,
                              wctparms *    parms,
                              int           ecount_same,
                              int           ecount_differ);
void iterate_dd(PricerSolver *solver);
void iterate_zdd(PricerSolver *solver);
size_t get_datasize(PricerSolver *solver);
void copy_solver(PricerSolver **dest, PricerSolver *src);
int add_one_conflict(
    PricerSolver *solver, wctparms *parms, int v1, int v2, int same);
int init_tables(PricerSolver *solver);
size_t get_numberrows_zdd(PricerSolver *solver);
size_t get_numberrows_bdd(PricerSolver *solver);
void set_release_due_time(PricerSolver *solver, Job *jobarray);

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
} data_status;

typedef struct wctdata wctdata;
struct wctdata {
    // The id and depth of the node in the B&B tree
    int id;
    int depth;

    data_status status;

    // The job information
    int njobs;
    int nmachines;
    // int *duetime;
    // int *releasetime;
    // int *duration;
    // int *weights;
    int *orig_node_ids;
    // data for meta heuristic
    Job *jobarray;
    int  H_max;
    int  H_min;

    // The column generation lp information
    wctlp * LP;
    double *x;
    double *coef;
    double *pi;
    // PricerSolver
    PricerSolver *solver;
    // Colorset(Assignments)
    int          ccount;
    Scheduleset *cclasses;
    int          dzcount;
    int          gallocated;
    Scheduleset *newsets;
    int          nnewsets;

    int  kpc_pi_scalef;
    int  kpc_pi_scalef_heur;
    int *kpc_pi;
    int *kpc_pi_heur;

    int     lower_bound;
    int     upper_bound;
    int     lower_scaled_bound;
    double  partial_sol;
    double  dbl_safe_lower_bound;
    double  dbl_est_lower_bound;
    double  dbl_est_lower_bound_heur;
    double  LP_lower_bound;
    double  LP_lower_bound_heur;
    double  LP_lower_bound_dual;
    double  LP_lower_bound_BB;
    double *rhs;
    int     nnonimprovements;
    /** Wentges smoothing technique */
    double *pi_in;
    double *pi_out;
    double *pi_sep;
    double *subgradient;
    double *subgradient_in;
    double  eta_in;
    double  eta_out;
    double  eta_sep;
    double  alpha;
    int     update;

    // Best Solution
    Scheduleset *bestcolors;
    int          besttotwct;
    int          nbbest;

    const Scheduleset *debugcolors;
    int                ndebugcolors;
    int                opt_track;

    // maxiterations and retireage
    int maxiterations;
    int retirementage;

    /** Branching strategies */
    int choose;
    /** conflict */
    int *    elist_same;
    int      ecount_same;
    int *    elist_differ;
    int      ecount_differ;
    wctdata *same_children;
    int      nsame;
    wctdata *diff_children;
    int      ndiff;
    int      v1, v2;
    /** ahv branching */
    wctdata *duetime_child;
    int      nduetime;
    wctdata *releasetime_child;
    int      nreleasetime;
    int      branch_job;
    int      completiontime;
    /** wide branching conflict */
    int *     v1_wide;
    int *     v2_wide;
    int       nb_wide;
    wctdata **same_children_wide;
    wctdata **diff_children_wide;

    wctdata *parent;

    char pname[MAX_PNAME_LEN];
};

/**
 * wct problem data type
 */

typedef enum {
    no_sol = 0,
    lp_feasible = 1,
    feasible = 2,
    meta_heur = 3,
    optimal = 4
} problem_status;

typedef struct wctproblem wctproblem;

struct wctproblem {
    wctparms parms;
    wctdata  root_pd;
    /** Job data */
    Job *jobarray;
    /** Summary of jobs */
    int njobs;
    int psum;
    int pmax;
    int pmin;
    int dmax;
    int dmin;
    int T;
    int off;
    /** EDD order of the jobs */
    Job **ojobarray;
    /** nmachines */
    int nmachines;

    int    nwctdata;
    int    global_upper_bound;
    int    first_upper_bound;
    int    global_lower_bound;
    int    first_lower_bound;
    double first_rel_error;
    double rel_error;
    int    maxdepth;

    problem_status status;

    /* All stable sets*/
    Scheduleset *initsets;
    int          nbinitsets;
    int          gallocated;
    /* Best Solution*/
    solution *   opt_sol;
    Scheduleset *bestschedule;
    int          nbestschedule;
    /*heap variables*/
    BinomialHeap *br_heap_a;
    GPtrArray *   unexplored_states;
    GQueue *      non_empty_level_pqs;
    unsigned int  last_explored;
    int           mult_key;
    int           found;
    int           nb_explored_nodes;
    int           nb_generated_col;
    int           nb_generated_col_root;
    /*Cpu time measurement*/
    CCutil_timer tot_build_dd;
    CCutil_timer tot_cputime;
    CCutil_timer tot_branch_and_bound;
    CCutil_timer tot_lb;
    CCutil_timer tot_lb_lp_root;
    CCutil_timer tot_lb_lp;
    CCutil_timer tot_pricing;
    CCutil_timer tot_scatter_search;
    double       real_time;
};

/*Initialization and free memory for the problem*/
void wctproblem_init(wctproblem *problem);
void wctproblem_free(wctproblem *problem);

/*Initialize pmc data*/
void wctdata_init(wctdata *pd);
int set_id_and_name(wctdata *pd, int id, const char *fname);
int wctdata_init_unique(wctdata *pd, int id, const char *name);
int lp_build(wctdata *pd);

/*Free the wctdata*/
void lpwctdata_free(wctdata *pd);
void children_data_free(wctdata *pd);
void wctdata_free(wctdata *pd);
void add_all_sets(PricerSolver *solver, wctdata *pd);
#ifdef __cplusplus
}
#endif

#endif
