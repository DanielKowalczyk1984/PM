// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "Parms.h"
#include <docopt/docopt.h>                       // for value, docopt_parse
#include <fmt/core.h>                            // for print
#include <array>                                 // for array
#include <cstddef>                               // for size_t
#include <fstream>                               // for ifstream
#include <nlohmann/json.hpp>                     // for json
#include <range/v3/iterator/basic_iterator.hpp>  // for operator-, operator!=
#include <range/v3/range/conversion.hpp>         // for to_container::fn
#include <range/v3/view/drop.hpp>                // for drop, drop_fn
#include <range/v3/view/transform.hpp>           // for transform_view, tran...
#include <regex>                                 // for regex_search, match_...
#include <span>                                  // for span
#include <string>                                // for allocator, string, stod
#include "Usage.hpp"                             // for USAGE
#include "orutils/util.h"                        // for dbg_lvl, program_header

const size_t TIME_LIMIT = 7200;
const double ALPHA_STAB_INIT = 0.8;

NLOHMANN_JSON_SERIALIZE_ENUM(PricingSolver,
                             {{bdd_solver_simple, "BddForward"},
                              {bdd_solver_cycle, "BddForwardCycle"},
                              {bdd_solver_backward_simple, "BddBackward"},
                              {bdd_solver_backward_cycle, "BddBackwardCycle"},
                              {zdd_solver_simple, "ZddForward"},
                              {zdd_solver_cycle, "ZddForwardCycle"},
                              {zdd_solver_backward_simple, "ZddBackward"},
                              {zdd_solver_backward_cycle, "ZddBackwardCycle"},
                              {dp_solver, "Time-Indexed"},
                              {ati_solver, "Arc-Time-Indexed"},
                              {dp_bdd_solver, "Hybrid"}})

NLOHMANN_JSON_SERIALIZE_ENUM(StabTechniques,
                             {
                                 {no_stab, "NoStabilization"},
                                 {stab_wentgnes, "WentgnesStab"},
                                 {stab_dynamic, "DynamicStab"},
                                 {stab_hybrid, "HybridStab"},
                             })

NLOHMANN_JSON_SERIALIZE_ENUM(BBExploreStrategy,
                             {{min_bb_explore_strategy, "dfs"},
                              {bb_dfs_strategy, "dfs"},
                              {bb_bfs_strategy, "bfs"},
                              {bb_brfs_strategy, "brfs"},
                              {bb_cbfs_strategy, "cbfs"}})

NLOHMANN_JSON_SERIALIZE_ENUM(Scoring_Parameter,
                             {{min_scoring_parameter, "ProductScoring"},
                              {product_scoring_parameter, "ProductScoring"},
                              {min_function_scoring_parameter, "MinFunction"},
                              {max_function_scoring_parameter, "MaxFunction"},
                              {weighted_sum_scoring_parameter, " WeightedSum"},
                              {weighted_product_scoring_parameter,
                               "WeightedProduct"}})

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
      use_heuristic(true),
      use_mip_solver(false),
      refine_bdd(false),
      enumerate(false),
      pruning_test(false),
      suboptimal_duals(false),
      reduce_cost_fixing(true),
      stab_technique(min_stab),
      print_csv(false),
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

    auto args = docopt::docopt(
        USAGE,
        tmp_char | ranges::views::drop(1) |
            ranges::views::transform(
                [](auto tmp) -> std::string { return std::string(tmp); }) |
            ranges::to_vector,
        true, "PM 0.1");

    if (!args["--json"].asBool()) {
        bb_node_limit = static_cast<int>(args["--node_limit"].asLong());
        branching_cpu_limit = static_cast<size_t>(args["--cpu_limit"].asLong());
        nb_iterations_rvnd = static_cast<int>(args["--nb_rvnb_it"].asLong());
        pricing_solver =
            static_cast<PricingSolver>(args["--pricing_solver"].asLong());
        strong_branching =
            static_cast<int>(args["--strong_branching"].asLong());

        bb_explore_strategy = static_cast<BBExploreStrategy>(
            static_cast<int>(args["--branching_strategy"].asLong()));
        stab_technique = static_cast<StabTechniques>(
            static_cast<int>(args["--stab_method"].asLong()));

        alpha = std::stod(args["--alpha"].asString());
        branching_point = std::stod(args["--branching_point"].asString());

        reduce_cost_fixing = !(args["--no_rc_fixing"].asBool());
        enumerate = args["--enumerate"].asBool();
        print_csv = args["--print_csv"].asBool();
        pruning_test = args["--pruning_test"].asBool();
        refine_bdd = args["--refinement"].asBool();
        suboptimal_duals = args["--suboptimal_duals"].asBool();
        use_heuristic = !(args["--no_heuristic"].asBool());
        use_mip_solver = args["--use_mip_solver"].asBool();
        parms_set_scoring_function(
            static_cast<int>(args["--scoring_function"].asLong()));

    } else {
        nlohmann::json j;
        std::ifstream  file(args["<json_file>"].asString());
        file >> j;
        j.at("bb_explore_strategy").get_to(bb_explore_strategy);
        j.at("scoring_parameter").get_to(scoring_parameter);
        j.at("scoring_value").get_to(scoring_value);
        j.at("strong_branching").get_to(strong_branching);
        j.at("bb_node_limit").get_to(bb_node_limit);
        j.at("nb_iterations_rvnd").get_to(nb_iterations_rvnd);
        j.at("branching_cpu_limit").get_to(branching_cpu_limit);
        j.at("alpha").get_to(alpha);
        j.at("branching_point").get_to(branching_point);
        j.at("pricing_solver").get_to(pricing_solver);
        j.at("use_heuristic").get_to(use_heuristic);
        j.at("use_mip_solver").get_to(use_mip_solver);
        j.at("refine_bdd").get_to(refine_bdd);
        j.at("enumerate").get_to(enumerate);
        j.at("pruning_test").get_to(pruning_test);
        j.at("suboptimal_duals").get_to(suboptimal_duals);
        j.at("reduce_cost_fixing").get_to(reduce_cost_fixing);
        j.at("stab_technique").get_to(stab_technique);
        j.at("print_csv").get_to(print_csv);
        parms_set_scoring_function(scoring_parameter);
    }

    /** Determine the name of the instance */
    auto file_name = args["FILE"].asString();
    jobfile = file_name;
    pname = find_match(file_name);
    /** Set the number of machines */
    nb_machines = static_cast<size_t>(args["NB"].asLong());

    return val;
}
