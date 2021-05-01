#include "Parms.h"
#include <docopt/docopt.h>
#include <fmt/core.h>
#include <algorithm>
#include <cstddef>
#include <regex>
#include <string>
#include "util.h"
#include "wctprivate.h"

const size_t TIME_LIMIT = 7200;
const double ALPHA_STAB_INIT = 0.8;

Parms::Parms()
    : init_upper_bound(INT_MAX),
      bb_explore_strategy(min_bb_explore_strategy),
      scoring_parameter(min_scoring_parameter),
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

Parms::Parms(int argc, const char** argv) : Parms() {
    program_header(argc, argv);

    parse_cmd(argc, argv);

    if (dbg_lvl() > 1) {
        fmt::print("Debugging turned on\n");
    }
}

int Parms::parms_set_file(std::string const& fname) {
    jobfile = std::string(fname);
    return 0;
}

int Parms::parms_set_pname(std::string const& fname) {
    pname = fname;
    return 0;
}

int Parms::parms_set_branching_cpu_limit(size_t limit) {
    branching_cpu_limit = limit;
    return 0;
}

int Parms::parms_set_alpha(double _alpha) {
    alpha = _alpha;
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

int Parms::parms_set_stab_technique(int _stab_technique) {
    stab_technique = static_cast<StabTechniques>(_stab_technique);
    return 0;
}

int Parms::parms_set_print(int _print) {
    print = _print;
    return 0;
}

int Parms::parms_set_pricing_solver(int solver) {
    pricing_solver = solver;
    return 0;
}

int Parms::parms_set_scoring_function(int scoring) {
    scoring_parameter = static_cast<Scoring_Parameter>(scoring);
    switch (scoring_parameter) {
        case min_scoring_parameter:
            scoring_function = [](double left, double right) {
                return std::max(left, EPS) * std::max(right, EPS);
            };

            break;
        case min_function_scoring_parameter:
            scoring_function = [](double left, double right) {
                return std::min(left, right);
            };
            break;
        case weighted_sum_scoring_parameter:
            scoring_function = [](double left, double right) {
                return (1.0 - mu) * std::max(left, right) +
                       mu * std::min(left, right);
            };
            break;
        case weighted_product_scoring_parameter:
            scoring_function = [](double left, double right) {
                return std::pow(std::max(left - 1.0, EPS), beta[0]) *
                       std::pow(std::max(right - 1.0, EPS), beta[1]);
            };
            break;
    }

    return 0;
}

static std::string find_match(std::string const& _instance_file) {
    std::regex  regexp{"^.*(wt[0-9]*_[0-9]*).*$"};
    std::smatch match{};
    std::regex_search(_instance_file, match, regexp);

    if (match.size() != 2) {
        return std::string("unknown_problem");
    } else {
        return match[1];
    }
}

int Parms::parse_cmd(int argc, const char** argv) {
    int val = 0;

    auto args =
        docopt::docopt_parse(USAGE, {argv + 1, argv + argc}, true, "PM 0.1");

    /** Set CPU limit for branching */
    parms_set_branching_cpu_limit(
        static_cast<size_t>(args["--cpu_limit"].asLong()));
    /** Set number of iterations in rvnd */
    parms_set_nb_iterations_rvnd(
        static_cast<int>(args["--nb_rvnb_it"].asLong()));
    /** Print statistics to csv files */
    parms_set_print(args["--print_csv"].asBool());
    /** Use MIP solver or not */
    parms_set_mip_solver(args["--mip_solver"].asBool());
    /** Use reduced cost fixing */
    parms_set_reduce_cost(!(args["--no_rc_fixing"].asBool()));
    /** Use heuristic or not */
    parms_set_use_heuristic(!(args["--no_heuristic"].asBool()));
    /** Set the pricing solver */
    parms_set_pricing_solver(
        static_cast<int>(args["--pricing_solver"].asLong()));
    /** Set the stabilization method */
    parms_set_stab_technique(static_cast<int>(args["--stab_method"].asLong()));
    parms_set_scoring_function(
        static_cast<int>(args["--scoring_function"].asLong()));
    parms_set_branchandbound(args["--no_branch_and_bound"].asBool());
    parms_set_strong_branching(!(args["--no_strong_branching"].asBool()));
    parms_set_alpha(std::stod(args["--alpha"].asString()));
    parms_set_bb_explore_strategy(
        static_cast<int>(args["--branching_strategy"].asLong()));
    parms_set_bb_node_limit(static_cast<int>(args["--node_limit"].asLong()));
    /** Determine the name of the instance */
    auto file_name = args["FILE"].asString();
    parms_set_file(file_name);
    parms_set_pname(find_match(file_name));
    /** Set the number of machines */
    parms_set_nb_machines(static_cast<int>(args["NB"].asLong()));

    return val;
}
