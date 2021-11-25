// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

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
#include <boost/chrono/duration.hpp>             // for duration::den
#include <boost/timer/timer.hpp>                 // for nano_second_types
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
#include "DebugLvl.hpp"                          // for debug_lvl
#include "Usage.hpp"                             // for USAGE
#include "orutils/util.h"                        // for program_header

const boost::timer::nanosecond_type TIME_LIMIT = 7200;
const double                        ALPHA_STAB_INIT = 0.8;

Parms::Parms()
    : bb_explore_strategy(min_bb_explore_strategy),
      scoring_parameter(min_scoring_parameter),
      scoring_value(min_scoring_value),
      strong_branching(),
      bb_node_limit(0),
      nb_iterations_rvnd(3),
      branching_cpu_limit(TIME_LIMIT * boost::chrono::nanoseconds::period::den),
      use_cpu_time{false},
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

    if (debug_lvl(1)) {
        fmt::print("Debugging turned on\n");
    }
}

void Parms::parms_set_scoring_function(int scoring) {
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

void to_json(nlohmann::json& j, const Parms& p) {
    j = nlohmann::json{{"bb_explore_strategy", p.bb_explore_strategy},
                       {"scoring_parameter", p.scoring_parameter},
                       {"scoring_value", p.scoring_value},
                       {"strong_branching", p.strong_branching},
                       {"bb_node_limit", p.bb_node_limit},
                       {"nb_iterations_rvnd", p.nb_iterations_rvnd},
                       {"branching_cpu_limit", p.branching_cpu_limit},
                       {"use_cpu_time", p.use_cpu_time},
                       {"alpha", p.alpha},
                       {"branching_point", p.branching_point},
                       {"pricing_solver", p.pricing_solver},
                       {"use_heuristic", p.use_heuristic},
                       {"use_mip_solver", p.use_mip_solver},
                       {"refine_bdd", p.refine_bdd},
                       {"enumerate", p.enumerate},
                       {"pruning_test", p.pruning_test},
                       {"suboptimal_duals", p.suboptimal_duals},
                       {"reduce_cost_fixing", p.reduce_cost_fixing},
                       {"stab_technique", p.stab_technique},
                       {"print_csv", p.print_csv}};
}

void from_json(const nlohmann::json& j, Parms& p) {
    j.at("bb_explore_strategy").get_to(p.bb_explore_strategy);
    j.at("scoring_parameter").get_to(p.scoring_parameter);
    j.at("scoring_value").get_to(p.scoring_value);
    j.at("strong_branching").get_to(p.strong_branching);
    j.at("bb_node_limit").get_to(p.bb_node_limit);
    j.at("nb_iterations_rvnd").get_to(p.nb_iterations_rvnd);
    j.at("branching_cpu_limit").get_to(p.branching_cpu_limit);
    p.branching_cpu_limit *= boost::chrono::nanoseconds::period::den;
    j.at("use_cpu_time").get_to(p.use_cpu_time);
    j.at("alpha").get_to(p.alpha);
    j.at("branching_point").get_to(p.branching_point);
    j.at("pricing_solver").get_to(p.pricing_solver);
    j.at("use_heuristic").get_to(p.use_heuristic);
    j.at("use_mip_solver").get_to(p.use_mip_solver);
    j.at("refine_bdd").get_to(p.refine_bdd);
    j.at("enumerate").get_to(p.enumerate);
    j.at("pruning_test").get_to(p.pruning_test);
    j.at("suboptimal_duals").get_to(p.suboptimal_duals);
    j.at("reduce_cost_fixing").get_to(p.reduce_cost_fixing);
    j.at("stab_technique").get_to(p.stab_technique);
    j.at("print_csv").get_to(p.print_csv);
}

void Parms::parse_cmd(int argc, const char** argv) {
    std::span tmp_char{argv, static_cast<size_t>(argc)};

    auto args = docopt::docopt(
        USAGE,
        tmp_char | ranges::views::drop(1) |
            ranges::views::transform(
                [](auto tmp) -> std::string { return std::string(tmp); }) |
            ranges::to_vector,
        true, "PM 0.1");

    if (!args["--json"].asBool()) {
        /**
         * @brief
         * All integer parameters
         */
        bb_node_limit = static_cast<int>(args["--node_limit"].asLong());
        branching_cpu_limit = static_cast<boost::timer::nanosecond_type>(
                                  args["--cpu_limit"].asLong()) *
                              boost::chrono::nanoseconds::period::den;
        nb_iterations_rvnd = static_cast<int>(args["--nb_rvnb_it"].asLong());
        strong_branching =
            static_cast<int>(args["--strong_branching"].asLong());

        /**
         * @brief
         * All the enumeration variables
         */
        pricing_solver =
            static_cast<PricingSolver>(args["--pricing_solver"].asLong());
        bb_explore_strategy = static_cast<BBExploreStrategy>(
            static_cast<int>(args["--branching_strategy"].asLong()));
        stab_technique = static_cast<StabTechniques>(
            static_cast<int>(args["--stab_method"].asLong()));

        /**
         * @brief
         * All the double parameters
         */
        alpha = std::stod(args["--alpha"].asString());
        branching_point = std::stod(args["--branching_point"].asString());

        /**
         * @brief all the bool parameters
         *
         */
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
        *this = j.get<Parms>();

        auto aux_dbg = -1;
        j.at("debug_lvl").get_to(aux_dbg);
        set_debug_lvl(aux_dbg);
        parms_set_scoring_function(scoring_parameter);
    }

    /** Determine the name of the instance */
    auto file_name = args["FILE"].asString();
    jobfile = file_name;
    pname = find_match(file_name);
    /** Set the number of machines */
    nb_machines = static_cast<size_t>(args["NB"].asLong());
}
