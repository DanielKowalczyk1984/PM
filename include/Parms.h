////////////////////////////////////////////////////////////////
//                                                            //
//  wctparms.h                                                //
//  wct                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef INCLUDE_WCTPARMS_H_
#define INCLUDE_WCTPARMS_H_
#include <cstddef>
#include <functional>
#include <string>

static const std::string USAGE =
    R"(PM.

Usage:
  bin/PM [-s <sn> -S <kn> -pmBRZHMd -n <nl> -b <br> -a <ln> -l <x> -f <y> -c <x> --alpha <mn> --branching_point <brp> --refinement] FILE NB
  bin/PM (-h | --help)
  bin/PM --version

Arguments:
  FILE  Path to the instance file
  NB    Number of machines

Options:
  -h --help                     Show this screen.
  --version                     Show version.
  -d --debug                    Turn on the debugging.
  -s --scoring_function=<sn>    Set scoring function branching[default: 0]
  -S --stab_method=<kn>         Stabilization technique: 0 = no stabilization, 1 = stabilization wentgnes, 2 = stabilization dynamic[default: 1].
  -a --pricing_solver=<ln>      Set pricing solver: 0 = bdd backward cycle, 1 = bdd forward simple, 2 = bdd forward cycle,
                                                    3 = bdd backward simple, 4 = bdd backward cycle, 5 = zdd forward simple,
                                                    6 = zdd forward cycle, 7 = zdd backward simple, 8 = zdd backward cycle,
                                                    9 = TI solver, 10 = arc-TI solver, 11 = hybrid model TI and bdd backward cycle[default: 4].
  -l --cpu_limit=<x>            Cpu time limit for branch and bound method[default: 7200].
  -f --nb_rvnb_it=<y>           Number of iterations in RVND[default: 5000].
  --alpha=<mn>                  Stabilization factor[default: 0.8].
  --branching_point=<brp>       Branching point[default: 0.2].
  -p --print_csv                Print csv-files.
  -r --refinement               Refine decision diagram.
  -m --mip_solver               Use mip solver to solve the original formulation.
  -R --no_rc_fixing             Don't apply reduce cost fixing.
  -H --no_heuristic             Don't apply heuristic.
  -B --no_branch_and_bound      Don't apply branch-and-bound.
  -c --strong_branching=<sb>    Don't apply strong branching[default: 20].
  -M --use_mip_solver           Use MIP solver.
  -b --branching_strategy=<br>  Set branch-and-bound exploration strategy: 0 = DFS, 1 = BFS, 2 = BrFS, 3 = CBFS[default: 0].
  -n --node_limit=<nl>          Set a limit on the number of nodes that can be explored.[default: 0]. Default meaning that all nodes should be explored.
)";

enum BBNodeSelection {
    min_search_strategy = 0,
    no_branching = min_search_strategy,
    min_lb_strategy = 1,
    dfs_strategy = 2,
    max_strategy = 3,
};

enum PricingSolver {
    bdd_solver_simple = 1,
    bdd_solver_cycle = 2,
    bdd_solver_backward_simple = 3,
    bdd_solver_backward_cycle = 4,
    zdd_solver_simple = 5,
    zdd_solver_cycle = 6,
    zdd_solver_backward_simple = 7,
    zdd_solver_backward_cycle = 8,
    dp_solver = 9,
    ati_solver = 10,
    dp_bdd_solver = 11,
};

enum BranchandBound {
    min_branch_and_bound = 0,
    no_branch_and_bound = 1,
    yes_branch_and_bound = min_branch_and_bound,
};

enum StabTechniques {
    no_stab = 0,
    stab_wentgnes = 1,
    stab_dynamic = 2,
    stab_hybrid = 3,
    min_stab = stab_wentgnes,
};

enum print {
    min_print_size = 0,
    use_print = 1,
};

enum BBBranchStrategy {
    min_bb_strategy = 0,
    conflict_strategy = min_bb_strategy,
    ahv_strategy = 1,
    cbfs_conflict_strategy = 2,
    cbfs_ahv_strategy = 3,
};

enum BBExploreStrategy {
    min_bb_explore_strategy = 0,
    bb_dfs_strategy = min_bb_explore_strategy,
    bb_bfs_strategy = 1,
    bb_brfs_strategy = 2,
    bb_cbfs_strategy = 3,
};

enum MIP_solver {
    min_mip_solver = 0,
    no_mip_solver = min_mip_solver,
    use_mip_solver = 1,
};

enum Scoring_Parameter {
    min_scoring_parameter = 0,
    product_scoring_parameter = min_scoring_parameter,
    min_function_scoring_parameter = 1,
    weighted_sum_scoring_parameter = 2,
    weighted_product_scoring_parameter = 3,
};

enum reduced_cost_fixing_param {
    min_reduced_cost = 1,
    yes_reduced_cost = min_reduced_cost,
    no_reduced_cost = 0,
};

enum use_heuristic {
    min_use_heuristic = 1,
    yes_use_heuristic = min_use_heuristic,
    no_use_heuristic = 0,
};

struct Parms {
    /**
     * General parameters
     */
    int                                   init_upper_bound;
    enum BBExploreStrategy                bb_explore_strategy;
    enum Scoring_Parameter                scoring_parameter;
    int                                   strong_branching;
    int                                   bb_node_limit;
    int                                   nb_iterations_rvnd;
    size_t                                branching_cpu_limit;
    double                                alpha;
    double                                branching_point;
    int                                   pricing_solver;
    int                                   mip_solver;
    int                                   use_heuristic;
    std::function<double(double, double)> scoring_function;
    bool                                  use_mip_solver;
    bool                                  refine_bdd;

    enum reduced_cost_fixing_param reduce_cost_fixing;

    /**
     * column generation
     */
    int                 branchandbound;
    enum StabTechniques stab_technique;
    int                 print;

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

    /*Functions for setting some parameters*/
    int parms_set_branching_cpu_limit(size_t limit);
    int parms_set_alpha(double alpha);
    int parms_set_strong_branching(int strong);
    int parms_set_mip_solver(int usage);
    int parms_set_use_heuristic(int usage);
    int parms_set_reduce_cost(int usage);
    int parms_set_nb_iterations_rvnd(int nb_sol);

    /**
     * Column Generation parameters
     */
    int parms_set_pricing_solver(int solver);
    int parms_set_stab_technique(int stab_technique);

    /**
     * Branch-and-Bound parameters
     */
    int parms_set_branchandbound(int bound);
    int parms_set_bb_explore_strategy(int strategy);
    int parms_set_bb_node_limit(int node_limit);
    int parms_set_print(int print);
    int parms_set_scoring_function(int scoring);

    /*Functions for defining the filesname*/
    int parms_set_file(std::string const& fname);
    int parms_set_pname(std::string const& fname);
    int parms_set_nb_machines(int nb_machines);
    int parms_set_use_mip_solver(bool _use_solver);
    /**
    @brief Parse the commands with docopt
     *
     * @param argc
     * @param argv
     * @param parms
     * @return int
     */
    int parse_cmd(int argc, const char** argv);

   private:
    static constexpr double                mu = 5.0 / 6.0;
    static constexpr std::array<double, 2> beta = {1.5, 0.5};
    static constexpr double                TargetBrTimeValue = 0.45;
};

#endif  // INCLUDE_WCTPARMS_H_
