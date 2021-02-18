#include "parms.h"
#include <docopt/docopt.h>
#include <limits.h>
#include <string.h>
#include <util.h>
#include <regex>
#include <string>
#include <vector>

const double TIME_LIMIT = 7200.0;
const double ALPHA_STAB_INIT = 0.8;

parms::parms()
    : init_upper_bound(INT_MAX),
      bb_explore_strategy(min_bb_explore_strategy),
      use_strong_branching(min_strong_branching),
      bb_node_limit(0),
      nb_iterations_rvnd(3),
      branching_cpu_limit(TIME_LIMIT),
      alpha(ALPHA_STAB_INIT),
      pricing_solver(bdd_solver_backward_cycle),
      mip_solver(min_mip_solver),
      use_heuristic(min_use_heuristic),
      reduce_cost_fixing(min_reduced_cost),
      branchandbound(min_branch_and_bound),
      stab_technique(min_stab),
      print(min_print_size),
      jobfile(),
      pname(),
      nb_jobs(0),
      nb_machines(0) {}

parms::parms(int argc, const char** argv) : parms() {
    parse_cmd(argc, argv);
}

int Parms::parms_set_file(std::string const& fname) {
    jobfile = std::string(fname);
    return 0;
}

int Parms::parms_set_pname(std::string const& fname) {
    pname = fname;
    return 0;
}

int Parms::parms_set_branching_cpu_limit(double limit) {
    branching_cpu_limit = limit;
    return 0;
}

int Parms::parms_set_branching_strategy(int strategy) {
    // bb_branch_strategy = strategy;
    return 0;
}

int Parms::parms_set_alpha(double alpha) {
    alpha = alpha;
    return 0;
}

int Parms::parms_set_mip_solver(int usage) {
    mip_solver = usage;
    return 0;
}

int Parms::parms_set_use_heuristic(int usage) {
    use_heuristic = usage;
    return 0;
}
int Parms::parms_set_reduce_cost(int usage) {
    reduce_cost_fixing = static_cast<reduced_cost_fixing_param>(usage);
    return 0;
}

int Parms::parms_set_strong_branching(int strong) {
    use_strong_branching = strong;
    return 0;
}

int Parms::parms_set_nb_machines(int _nb_machines) {
    nb_machines = _nb_machines;
    return 0;
}

int Parms::parms_set_nb_iterations_rvnd(int nb_iterations) {
    nb_iterations_rvnd = nb_iterations;
    return 0;
}

int Parms::parms_set_branchandbound(int bound) {
    branchandbound = bound;
    return 0;
}

int Parms::parms_set_bb_explore_strategy(int strategy) {
    bb_explore_strategy = static_cast<BBExploreStrategy>(strategy);
    return 0;
}

int Parms::parms_set_bb_node_limit(int node_limit) {
    bb_node_limit = node_limit;
    return 0;
}

int Parms::parms_set_stab_technique(int stab_technique) {
    stab_technique = static_cast<stab_techniques>(stab_technique);
    return 0;
}

int Parms::parms_set_print(int print) {
    print = print;
    return 0;
}

int Parms::parms_set_pricing_solver(int solver) {
    pricing_solver = solver;
    return 0;
}

static const std::string USAGE =
    R"(PM.

Usage:
  bin/PM [-S <kn> -pmBRZH -n <nl> -b <br> -a <ln> -l <x> -f <y> -d --no_strong_branching --alpha <mn>] FILE NB
  bin/PM (-h | --help)
  bin/PM --version

Arguments:
  FILE  Path to the instance file
  NB    Number of machines

Options:
  -h --help                     Show this screen.
  --version                     Show version.
  -d --debug                    Turn on the debugging.
  -S --stab_method=<kn>         Stabilization technique: 0 = no stabilization, 1 = stabilization wentgnes, 2 = stabilization dynamic[default: 1].
  -a --pricing_solver=<ln>      Set pricing solver: 0 = bdd backward cycle, 1 = bdd forward simple, 2 = bdd forward cycle,
                                                    3 = bdd backward simple, 4 = bdd backward cycle, 5 = zdd forward simple,
                                                    6 = zdd forward cycle, 7 = zdd backward simple, 8 = zdd backward cycle,
                                                    9 = TI solver, 10 = arc-TI solver, 11 = hybrid model TI and bdd backward cycle[default: 4].
  -l --cpu_limit=<x>            Cpu time limit for branch and bound method[default: 7200].
  -f --nb_rvnb_it=<y>           Number of iterations in RVND[default: 4].
  --alpha=<mn>                  Stabilization factor[default: 0.8].
  -p --print_csv                Print csv-files.
  -m --mip_solver               Use mip solver to solve the original formulation.
  -R --no_rc_fixing             Don't apply reduce cost fixing.
  -H --no_heuristic             Don't apply heuristic.
  -B --no_branch_and_bound      Don't apply branch-and-bound.
  --no_strong_branching         Don't apply strong branching.
  -b --branching_strategy=<br>  Set branch-and-bound exploration strategy: 0 = DFS, 1 = BFS, 2 = BrFS, 3 = CBFS[default: 0].
  -n --node_limit=<nl>          Set a limit on the number of nodes that can be explored.[default: 0]. Default meaning that all nodes should be explored.
)";

static std::string find_match(std::string const& _instance_file) {
    std::regex  regexp{"^.*(wt[0-9]*_[0-9]*).*$"};
    std::smatch match{};
    // std::string s{_instance_file};
    std::regex_search(_instance_file, match, regexp);

    if (match.size() != 2) {
        return std::string("unknown_problem");
    } else {
        return match[1];
    }
}

int Parms::parse_cmd(int argc, const char** argv) {
    int val = 0;

    auto args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, "PM 0.1");

    /** Set CPU limit for branching */
    val = parms_set_branching_cpu_limit(
        static_cast<double>(args["--cpu_limit"].asLong()));
    /** Set number of iterations in rvnd */
    val = parms_set_nb_iterations_rvnd(
        static_cast<int>(args["--nb_rvnb_it"].asLong()));
    /** Print statistics to csv files */
    val = parms_set_print(args["--print_csv"].asBool());
    /** Use MIP solver or not */
    val = parms_set_mip_solver(args["--mip_solver"].asBool());
    /** Use reduced cost fixing */
    val = parms_set_reduce_cost(!(args["--no_rc_fixing"].asBool()));
    /** Use heuristic or not */
    val = parms_set_use_heuristic(!(args["--no_heuristic"].asBool()));
    /** Set the pricing solver */
    val = parms_set_pricing_solver(
        static_cast<int>(args["--pricing_solver"].asLong()));
    /** Set the stabilization method */
    val = parms_set_stab_technique(
        static_cast<int>(args["--stab_method"].asLong()));
    val = parms_set_branchandbound(args["--no_branch_and_bound"].asBool());
    val = parms_set_strong_branching(!(args["--no_strong_branching"].asBool()));
    val = parms_set_alpha(std::stod(args["--alpha"].asString()));
    val = parms_set_bb_explore_strategy(
        static_cast<int>(args["--branching_strategy"].asLong()));
    val = parms_set_bb_node_limit(
        static_cast<int>(args["--node_limit"].asLong()));
    /** Determine the name of the instance */
    auto file_name = args["FILE"].asString();
    val = parms_set_file(file_name);
    val = parms_set_pname(find_match(file_name));
    /** Set the number of machines */
    val = parms_set_nb_machines(static_cast<int>(args["NB"].asLong()));

    return val;
}
