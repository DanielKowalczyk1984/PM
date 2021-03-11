#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H

#include <bits/c++config.h>
#include <memory>
#include <vector>
#include "BranchBoundTree.hpp"
#include "Instance.h"
#include "MIP_defs.hpp"
#include "OptimalSolution.hpp"
#include "PricingStabilization.hpp"
#include "Solution.hpp"
#include "Statistics.h"
#include "binomial-heap.h"
// #include "interval.h"
#include "lp.h"
/**
 * problem data
 */
struct ScheduleSet;

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

constexpr int    NB_CUTS = 2000;
constexpr int    NB_CG_ITERATIONS = 1000000;
constexpr int    CLEANUP_ITERATION = 30;
constexpr double EPS = 1e-6;
constexpr double EPS_BOUND = 1e-9;

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

struct Problem {
    /** Different Parameters */
    Parms parms;
    /*Cpu time measurement + Statistics*/
    Statistics stat;
    /** Instance data*/
    Instance instance;

    std::unique_ptr<BranchBoundTree> tree;
    std::unique_ptr<NodeData>        root_pd;

    /** Summary of jobs */
    int nb_jobs;
    int nb_machines;

    int    global_upper_bound;
    int    global_lower_bound;
    double rel_error;
    int    root_upper_bound;
    int    root_lower_bound;
    double root_rel_error;

    problem_status status;

    /* Best Solution*/
    Sol opt_sol;

    /** All methods of problem class */
    int  to_screen();
    void to_csv();
    void solve();
    /** Heuristic related */
    void heuristic();
    /** Constructors */
    Problem(int argc, const char** argv);
    Problem(const Problem&) = delete;
    Problem(Problem&&) = delete;
    Problem& operator=(const Problem&) = delete;
    Problem& operator=(Problem&&) = delete;
    ~Problem();

    class ProblemException : public std::exception {
       public:
        ProblemException(const char* const msg = 0) { errmsg = msg; }

        const char* what(void) const throw() override { return (errmsg); }

       private:
        const char* errmsg;
    };
};

struct NodeData {
    // The id and depth of the node in the B&B tree
    int id;
    int depth;

    NodeDataStatus status;

    // The instance information
    const Instance& instance;
    int             nb_jobs;
    int             nb_machines;

    // The column generation lp information
    wctlp*              RMP;
    std::vector<double> lambda;

    std::vector<double> pi;
    std::vector<double> slack;
    std::vector<double> rhs;
    std::vector<double> lhs_coeff;
    std::vector<int>    id_row;
    std::vector<double> coeff_row;

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
    std::unique_ptr<PricerSolverBase> solver;

    // Columns
    int                                       zero_count;
    std::shared_ptr<ScheduleSet>              newsets;
    int                                       nb_new_sets;
    std::vector<int>                          column_status;
    std::vector<std::shared_ptr<ScheduleSet>> localColPool;

    int lower_bound;
    int upper_bound;

    double LP_lower_bound;
    double LP_lower_bound_dual;
    double LP_lower_bound_BB;
    double LP_lower_min;

    int nb_non_improvements;
    int iterations;

    /** Wentges smoothing technique */
    std::unique_ptr<PricingStabilizationBase> solver_stab;

    // Best Solution
    std::vector<std::shared_ptr<ScheduleSet>> best_schedule;
    int                                       best_objective;

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
    Sol*        opt_sol;
    std::string pname;

    explicit NodeData(Problem* problem);
    explicit NodeData(const Instance&);
    NodeData(const NodeData&) = delete;
    NodeData& operator=(const NodeData&) = delete;
    NodeData& operator=(NodeData&&) = delete;
    NodeData(NodeData&&) = default;
    ~NodeData();

    // void init_null();
    void lp_node_data_free();
    void temporary_data_free();
    void prune_duplicated_sets();
    // void add_solution_to_colpool(Solution*);
    void add_solution_to_colpool(const Sol&);
    // void add_solution_to_colpool_and_lp(Solution*);
    void add_solution_to_colpool_and_lp(const Sol&);

    int build_rmp();
    /** lowerbound.cpp */
    int  delete_unused_rows();
    int  delete_old_schedules();
    int  delete_infeasible_schedules();
    void make_pi_feasible_farkas_pricing();
    int  compute_objective();
    int  solve_relaxation();
    int  compute_lower_bound();
    int  check_schedules();
    int  print_x();

    /** PricerSolverWrappers.cpp */
    void build_solve_mip();
    void construct_lp_sol_from_rmp();
    void generate_cuts();
    int  delete_unused_rows_range(int first, int last);
    int  call_update_rows_coeff();
    bool check_schedule_set(ScheduleSet* set);
    void make_schedule_set_feasible(ScheduleSet* set);
    void get_mip_statistics(enum MIP_Attr c);

    /** StabilizationWrappers.cpp */
    int  solve_pricing();
    void solve_farkas_dbl();

    std::unique_ptr<NodeData> clone() const;

    int add_scheduleset_to_rmp(ScheduleSet* set);

   private:
    int add_lhs_scheduleset_to_rmp(ScheduleSet* set);

    /** lowerbound.cpp */
    int  grow_ages();
    void print_ages();

    /** StabilizationWrappers.cpp */
    template <typename T>
    int construct_sol(OptimalSolution<T>* sol);
};

#endif
