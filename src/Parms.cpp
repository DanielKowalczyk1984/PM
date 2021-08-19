#include "Parms.h"
#include <docopt/docopt.h>                       // for value, docopt_parse
#include <fmt/core.h>                            // for print
#include <algorithm>                             // for max, min
#include <array>                                 // for array
#include <climits>                               // for INT_MAX
#include <cmath>                                 // for pow
#include <cstddef>                               // for size_t
#include <map>                                   // for map
#include <range/v3/iterator/basic_iterator.hpp>  // for operator-, operator!=
#include <range/v3/range/conversion.hpp>         // for to_container::fn
#include <range/v3/view/drop.hpp>                // for drop, drop_fn
#include <range/v3/view/subrange.hpp>            // for subrange
#include <range/v3/view/transform.hpp>           // for transform_view, tran...
#include <range/v3/view/view.hpp>                // for operator|, view_closure
#include <regex>                                 // for regex_search, match_...
#include <span>                                  // for span
#include <string>                                // for allocator, string, stod
#include <vector>                                // for vector
#include "util.h"                                // for dbg_lvl, program_header

const size_t TIME_LIMIT = 7200;
const double ALPHA_STAB_INIT = 0.8;

Parms::Parms()
    : bb_explore_strategy(min_bb_explore_strategy),
      scoring_parameter(min_scoring_parameter),
      scoring_value(min_scoring_value),
      strong_branching(),
      bb_node_limit(0),
      nb_iterations_rvnd(3),
      branching_cpu_limit(TIME_LIMIT),
      alpha(ALPHA_STAB_INIT),
      branching_point(TargetBrTimeValue),
      pricing_solver(bdd_solver_backward_cycle),
      use_heuristic(min_use_heuristic),
      use_mip_solver(false),
      refine_bdd(false),
      enumerate(false),
      pruning_test(false),
      suboptimal_duals(false),
      reduce_cost_fixing(min_reduced_cost),
      stab_technique(min_stab),
      print_csv(min_print_size),
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

int Parms::parms_set_scoring_function(int scoring) {
    scoring_parameter = static_cast<Scoring_Parameter>(scoring);
    switch (scoring_parameter) {
        case min_scoring_parameter:
            scoring_function = [](const std::array<double, 2>& a) {
                return std::max(a[0], EPS) * std::max(a[1], EPS);
            };

            break;
        case min_function_scoring_parameter:
            scoring_function = [](const std::array<double, 2>& a) {
                return ranges::min(a);
            };
            break;
        case weighted_sum_scoring_parameter:
            scoring_function = [](const std::array<double, 2>& a) {
                return (1.0 - mu) * ranges::max(a) + mu * ranges::min(a);
            };
            break;
        case weighted_product_scoring_parameter:
            scoring_function = [](const std::array<double, 2>& a) {
                return std::pow(std::max(a[0] - 1.0, EPS), beta[0]) *
                       std::pow(std::max(a[1] - 1.0, EPS), beta[1]);
            };
            break;
        case max_function_scoring_parameter:
            scoring_function = [](const std::array<double, 2>& a) {
                return ranges::max(a);
            };
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

    std::span tmp_char{argv, static_cast<size_t>(argc)};

    auto args = docopt::docopt_parse(
        USAGE,
        tmp_char | ranges::views::drop(1) |
            ranges::views::transform(
                [](auto tmp) -> std::string { return std::string(tmp); }) |
            ranges::to_vector,
        true, "PM 0.1");

    bb_node_limit = static_cast<int>(args["--node_limit"].asLong());
    branching_cpu_limit = static_cast<size_t>(args["--cpu_limit"].asLong());
    nb_iterations_rvnd = static_cast<int>(args["--nb_rvnb_it"].asLong());
    pricing_solver = static_cast<int>(args["--pricing_solver"].asLong());
    strong_branching = static_cast<int>(args["--strong_branching"].asLong());

    bb_explore_strategy = static_cast<BBExploreStrategy>(
        static_cast<int>(args["--branching_strategy"].asLong()));
    stab_technique = static_cast<StabTechniques>(
        static_cast<int>(args["--stab_method"].asLong()));

    alpha = std::stod(args["--alpha"].asString());
    branching_point = std::stod(args["--branching_point"].asString());

    reduce_cost_fixing =
        static_cast<ReducedCostFixingParam>(!(args["--no_rc_fixing"].asBool()));
    enumerate = args["--enumerate"].asBool();
    print_csv = args["--print_csv"].asBool();
    pruning_test = args["--pruning_test"].asBool();
    refine_bdd = args["--refinement"].asBool();
    suboptimal_duals = args["--suboptimal_duals"].asBool();
    use_heuristic = !(args["--no_heuristic"].asBool());
    use_mip_solver = args["--use_mip_solver"].asBool();

    parms_set_scoring_function(
        static_cast<int>(args["--scoring_function"].asLong()));

    /** Determine the name of the instance */
    auto file_name = args["FILE"].asString();
    jobfile = std::string(file_name);
    pname = find_match(file_name);
    /** Set the number of machines */
    nb_machines = static_cast<int>(args["NB"].asLong());

    return val;
}
