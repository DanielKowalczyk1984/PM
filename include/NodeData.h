#ifndef __NODEDATA_H__
#define __NODEDATA_H__

#include <cstddef>

// #include "Instance.h"    // for Instance
// #include "Parms.h"       // for Parms
#include "Solution.hpp"  // for Sol
#include "Statistics.h"  // for Statistics
#include "lp.h"

/** Forward Declarations */

class PricingStabilizationBase;
class Problem;
struct PricerSolverBase;
struct Instance;
struct Parms;
struct Column;  // lines 18-18

struct NodeData {
    enum NodeDataStatus {
        initialized = 0,
        LP_bound_estimated = 1,
        LP_bound_computed = 2,
        submitted_for_branching = 3,
        infeasible = 4,
        finished = 5,
    };

    size_t depth;

    NodeDataStatus status;

    // The instance information
    const Parms&    parms;
    const Instance& instance;
    Statistics&     stat;
    Sol&            opt_sol;
    std::string     pname;

    size_t nb_jobs;
    size_t nb_machines;

    // The column generation lp information
    std::unique_ptr<wctlp, std::function<void(wctlp*)>> RMP;

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
    size_t max_nb_cuts;
    int    id_convex_constraint;
    int    id_assignment_constraint;
    int    id_valid_cuts;

    size_t id_art_var_convex;
    int    id_art_var_assignment;
    size_t id_art_var_cuts;
    size_t id_next_var_cuts;
    int    id_pseudo_schedules;

    // PricerSolver
    std::unique_ptr<PricerSolverBase> solver;

    // Columns
    int                                  zero_count;
    std::shared_ptr<Column>              newsets;
    int                                  nb_new_sets;
    std::vector<int>                     column_status;
    std::vector<std::shared_ptr<Column>> localColPool;

    int lower_bound;
    int upper_bound;

    double LP_lower_bound;
    double LP_lower_bound_dual;
    double LP_lower_bound_BB;
    double LP_lower_min;

    int    nb_non_improvements;
    size_t iterations;

    /** Wentges smoothing technique */
    std::unique_ptr<PricingStabilizationBase> solver_stab;

    // maxiterations and retireage
    int retirementage;

    /** Branching strategies */
    size_t branch_job;
    int    completiontime;
    int    less;

    /**
     * ptr to the parent node
     */
    explicit NodeData(Problem* problem);
    NodeData(const NodeData&);
    NodeData& operator=(const NodeData&) = delete;
    NodeData& operator=(NodeData&&) = delete;
    NodeData(NodeData&&) = delete;
    ~NodeData();

    void prune_duplicated_sets();
    void add_solution_to_colpool(const Sol&);

    int build_rmp();
    /** lowerbound.cpp */
    int  delete_unused_rows();
    int  delete_old_columns();
    int  delete_infeasible_columns();
    int  compute_objective();
    int  solve_relaxation();
    int  compute_lower_bound();
    int  estimate_lower_bound(size_t _iter);
    bool refinement();
    void make_pi_feasible_farkas_pricing();

    /** PricerSolverWrappers.cpp */
    void build_solve_mip();
    void construct_lp_sol_from_rmp();
    void generate_cuts();
    int  delete_unused_rows_range(int first, int last);
    int  call_update_rows_coeff();

    /** small getters */
    double get_score_value();

    /** StabilizationWrappers.cpp */
    int  solve_pricing();
    void solve_farkas_dbl();

    [[nodiscard]] std::unique_ptr<NodeData>  clone() const;
    [[nodiscard]] std::unique_ptr<NodeData>  clone(size_t _j,
                                                   int    _t,
                                                   bool   _left) const;
    std::array<std::unique_ptr<NodeData>, 2> create_child_nodes(size_t _j,
                                                                int    _t);

    std::array<std::unique_ptr<NodeData>, 2> create_child_nodes(size_t _j,
                                                                long   _t);
    int add_scheduleset_to_rmp(Column* set);

   private:
    int add_lhs_scheduleset_to_rmp(Column* set);

    /** lowerbound.cpp */
    int  grow_ages();
    void print_ages();

    static constexpr auto CLEANUP_ITERATION = 30;
    static constexpr auto EPS = 1e-6;
    static constexpr auto EPS_BOUND = 1e-9;
    static constexpr auto min_nb_del_row_ratio = 0.9;
    static constexpr auto NB_CG_ITERATIONS = 1000000UL;
    static constexpr auto NB_CUTS = 2000;
    static constexpr auto NB_ESTIMATE_IT = 20;
};

#endif  // __NODEDATA_H__