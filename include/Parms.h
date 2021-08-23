#ifndef __PARMS_H__
#define __PARMS_H__

#include <array>       // for array
#include <cstddef>     // for size_t
#include <functional>  // for function
#include <string>      // for allocator, string

enum BBNodeSelection {
    min_search_strategy = 0,
    no_branching = min_search_strategy,
    min_lb_strategy = 1,
    dfs_strategy = 2,
    max_strategy = 3,
};

enum PricingSolver {
    bdd_solver_simple = 0,
    bdd_solver_cycle = 1,
    bdd_solver_backward_simple = 2,
    bdd_solver_backward_cycle = 3,
    zdd_solver_simple = 4,
    zdd_solver_cycle = 5,
    zdd_solver_backward_simple = 6,
    zdd_solver_backward_cycle = 7,
    dp_solver = 8,
    ati_solver = 9,
    dp_bdd_solver = 10
};

enum StabTechniques {
    no_stab = 0,
    stab_wentgnes = 1,
    stab_dynamic = 2,
    stab_hybrid = 3,
    min_stab = stab_wentgnes,
};

enum BBExploreStrategy {
    min_bb_explore_strategy = 0,
    bb_dfs_strategy = min_bb_explore_strategy,
    bb_bfs_strategy = 1,
    bb_brfs_strategy = 2,
    bb_cbfs_strategy = 3,
};

enum Scoring_Parameter {
    min_scoring_parameter = 0,
    product_scoring_parameter = min_scoring_parameter,
    min_function_scoring_parameter = 1,
    weighted_sum_scoring_parameter = 2,
    weighted_product_scoring_parameter = 3,
    max_function_scoring_parameter = 4,
};

enum Scoring_Value {
    min_scoring_value = 0,
    lb_scoring_value = min_scoring_value,
    size_scoring_value = 1,
    nb_paths_scoring_value = 2,
};

struct Parms {
    enum BBExploreStrategy                              bb_explore_strategy;
    enum Scoring_Parameter                              scoring_parameter;
    enum Scoring_Value                                  scoring_value;
    int                                                 strong_branching;
    int                                                 bb_node_limit;
    int                                                 nb_iterations_rvnd;
    size_t                                              branching_cpu_limit;
    double                                              alpha;
    double                                              branching_point;
    enum PricingSolver                                  pricing_solver;
    int                                                 use_heuristic;
    std::function<double(const std::array<double, 2>&)> scoring_function;
    bool                                                use_mip_solver;
    bool                                                refine_bdd;
    bool                                                enumerate;
    bool                                                pruning_test;
    bool                                                suboptimal_duals;
    bool                                                reduce_cost_fixing;
    enum StabTechniques                                 stab_technique;

    /**
     * column generation
     */
    int print_csv;

    std::string jobfile;
    std::string pname;

    int nb_jobs;
    int nb_machines;
    Parms();
    Parms(int argc, const char** argv);

    Parms(const Parms&) = default;
    Parms(Parms&&) = default;
    Parms& operator=(const Parms&) = default;
    Parms& operator=(Parms&&) = default;
    ~Parms() = default;

    int parms_set_scoring_function(int scoring);
    int parse_cmd(int argc, const char** argv);

   private:
    static constexpr double                mu = 5.0 / 6.0;
    static constexpr std::array<double, 2> beta = {1.5, 0.5};
    static constexpr double                TargetBrTimeValue = 0.45;
    static constexpr auto                  EPS = 1e-6;
};

#endif  // __PARMS_H__