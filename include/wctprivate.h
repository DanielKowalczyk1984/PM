#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H

#include <binomial-heap.h>
#include <lp.h>
#include <util.h>
#include <solver.h>
#include <interval.h>

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
} data_status;

/**
 * problem data
 */
typedef struct wctproblem wctproblem;
/**
 * node data
 */
typedef struct wctdata wctdata;

struct wctdata {
    // The id and depth of the node in the B&B tree
    int id;
    int depth;

    data_status status;

    // The instance information
    int njobs;
    int nmachines;
    int *orig_node_ids;
    // data for meta heuristic
    GPtrArray *jobarray;
    int  H_max;
    int  H_min;
    /** data about the intervals */
    GPtrArray *local_intervals;
    GPtrArray *ordered_jobs;
    int **sump;

    // The column generation lp information
    wctlp *LP;
    wctlp *MIP;
    double *x;
    double *x_e;
    double *coef;
    double *pi;
    // PricerSolver
    PricerSolver *solver;

    // Columns
    // int          ccount;
    // scheduleset *cclasses;
    int          dzcount;
    // int          gallocated;
    scheduleset *newsets;
    int          nnewsets;
    int *cstat;
    GPtrArray *localColPool;

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
    int iterations;
    /** Wentges smoothing technique */
    double *pi_in;
    double dualdiffnorm;
    double *subgradient;
    double hybridfactor;
    double subgradientnorm;
    double  alpha;
    double alphabar;
    double beta;
    int k;
    int node_stab;
    int     hasstabcenter;
    double  eta_in;
    int inmispricingschedule;
    double subgradientproduct;
    double *pi_out;
    double *pi_sep;
    double *subgradient_in;
    double  eta_out;
    double  eta_sep;
    double reduced_cost;
    int     update;

    // Best Solution
    scheduleset *bestcolors;
    int          besttotwct;
    int          nbbest;

    const scheduleset *debugcolors;
    int                ndebugcolors;
    int                opt_track;

    // maxiterations and retireage
    int maxiterations;
    int retirementage;

    /** Branching strategies */
    int choose;
    /** conflict */
    int     *elist_same;
    int      ecount_same;
    int     *elist_differ;
    int      ecount_differ;
    wctdata *same_children;
    int      nsame;
    wctdata *diff_children;
    int      ndiff;
    Job      *v1, *v2;
    /** ahv branching */
    wctdata *duetime_child;
    int      nduetime;
    wctdata *releasetime_child;
    int      nreleasetime;
    int      branch_job;
    int      completiontime;
    /** wide branching conflict */
    int      *v1_wide;
    int      *v2_wide;
    int       nb_wide;
    wctdata **same_children_wide;
    wctdata **diff_children_wide;


    /**
     * ptr to the parent node
     */
    wctdata *parent;

    /**
     * ptr to the data overview
     */
    wctproblem *problem;

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


struct wctproblem {
    wctparms parms;
    wctdata  root_pd;
    /** Job data in EDD order */
    GPtrArray *g_job_array;
    /** Summary of jobs */
    int njobs;
    int psum;
    int pmax;
    int pmin;
    int dmax;
    int dmin;
    int H_min;
    int H_max;
    int off;
    /** nmachines */
    int nmachines;

    int    nwctdata;
    int    global_upper_bound;
    int    global_lower_bound;
    double rel_error;
    int    root_upper_bound;
    int    root_lower_bound;
    double root_rel_error;
    int    maxdepth;

    problem_status status;

    /* All partial schedules*/
    scheduleset *initsets;
    int          nbinitsets;
    int          gallocated;
    GPtrArray *ColPool;
    /** Maximum number of artificial columns */
    int maxArtificials;
    /** Actual number of artificial columns */
    int nArtificials;
    /* Best Solution*/
    solution    *opt_sol;
    scheduleset *bestschedule;
    int          nbestschedule;
    /*heap variables*/
    BinomialHeap *br_heap_a;
    GPtrArray    *unexplored_states;
    GQueue       *non_empty_level_pqs;
    unsigned int  last_explored;
    int           mult_key;
    int           found;
    int           nb_explored_nodes;
    int           nb_generated_nodes;
    int           nb_generated_col;
    int           nb_generated_col_root;
    /*Cpu time measurement*/
    CCutil_timer tot_build_dd;
    CCutil_timer tot_cputime;
    CCutil_timer tot_branch_and_bound;
    CCutil_timer tot_strong_branching;
    CCutil_timer tot_lb_root;
    CCutil_timer tot_lb;
    CCutil_timer tot_solve_lp;
    CCutil_timer tot_pricing;
    CCutil_timer tot_heuristic;

    double real_time_build_dd;
    double real_time_total;
    double real_time_branch_and_bound;
    double real_time_strong_branching;
    double real_time_lb_root;
    double real_time_lb;
    double real_time_solve_lp;
    double real_time_pricing;
    double real_time_heuristic;
    //double       real_time;
};

/*Initialization and free memory for the problem*/
void wctproblem_init(wctproblem *problem);
void wctproblem_free(wctproblem *problem);

/*Initialize pmc data*/
void wctdata_init(wctdata *pd, wctproblem *prob);
int set_id_and_name(wctdata *pd, int id, const char *fname);
int wctdata_init_unique(wctdata *pd, int id, const char *name);
int lp_build(wctdata *pd);

/*Free the wctdata*/
void lpwctdata_free(wctdata *pd);
void children_data_free(wctdata *pd);
void wctdata_free(wctdata *pd);
void temporary_data_free(wctdata *pd);

/**
 * solver zdd
 */
int evaluate_nodes(wctdata *pd);
int calculate_new_ordered_jobs(wctdata *pd);
int build_solve_mip(wctdata *pd);


/**
 * pricing algorithms
 */
int solve_pricing(wctdata *pd, wctparms *parms, int evaluate);
int solve_stab(wctdata *pd, wctparms *parms);
int solve_stab_dynamic(wctdata *pd, wctparms *parms);
int solve_stab_hybrid(wctdata *pd, wctparms *parms);
int solve_farkas_dbl(wctdata *pd);
int solve_farkas_dbl_DP(wctdata *pd);
#ifdef __cplusplus
}
#endif

#endif
