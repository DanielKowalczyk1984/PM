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
#include <fmt/chrono.h>         // for print
#include <cstdio>               // for fopen, fclose, FILE
#include <ctime>                // for localtime, time, time_t
#include <filesystem>           // for path, directory_entry
#include <functional>           // for function
#include <memory>               // for operator==, unique_ptr
#include <string>               // for string, basic_string
#include "BranchBoundTree.hpp"  // for BranchBoundTree
#include "Instance.h"           // for Instance
#include "Parms.h"              // for Parms
#include "Problem.h"            // for Problem
#include "Statistics.h"         // for Statistics, Statistics::bb_timer, Sta....

/**
 * @brief Print the results to a csv file with fmt.
 *
 */
void Problem::to_csv() {
    using ptr_file = std::unique_ptr<std::FILE, std::function<int(FILE*)>>;

    ptr_file file{};
    auto     result = std::time(nullptr);

    auto file_name =
        fmt::format("CG_overall_{:%Y_%m_%d}.csv", fmt::localtime(result));
    stat.suspend_timer(Statistics::cputime_timer);
    stat.real_time_total = getRealTime() - stat.real_time_total;
    auto path_file = std::filesystem::current_path() / file_name;
    std::filesystem::directory_entry tmp_entry_file{path_file};

    if (tmp_entry_file.exists()) {
        file = ptr_file(std::fopen(file_name.c_str(), "a"), &fclose);
    } else {
        file = ptr_file(std::fopen(file_name.c_str(), "w"), &fclose);
        if (file == nullptr) {
            throw ProblemException(
                fmt::format("We could not open file {} in {} at {}", file_name,
                            __FILE__, __LINE__)
                    .c_str());
        }
        fmt::print(
            file.get(),
            R"({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}
)",
            "n", "m", "NameInstance", "tot_real_time", stat.time_total.name(),
            stat.time_branch_and_bound.name(), stat.time_lb.name(),
            stat.time_lb_root.name(), stat.time_heuristic.name(),
            stat.time_build_dd.name(), stat.time_pricing.name(),
            stat.time_rc_fixing.name(), "rel_error", "global_lower_bound",
            "global_upper_bound", "first_rel_error", "global_upper_bound_root",
            "global_lowerbound_root", "nb_generated_col",
            "nb_generated_col_root", "first_size_graph",
            "size_after_reduced_cost", "mip_nb_vars", "mip_nb_constr",
            "mip_obj_bound", "mip_obj_bound_lp", "mip_rel_gap", "mip_run_time",
            "mip_status", "mip_nb_iter_simplex", "mip_nb_nodes",
            "nb_nodes_explored", "date", "nb_iterations_rvnd", "stabilization",
            "alpha", "pricing_solver", "strong_branching", "branching_point",
            "refinement", "pruning_test", "suboptimal_duals",
            "scoring_parameter", "scoring_value", "use_cpu_time");
    }

    fmt::print(file.get(),
               R"({},{},{},{:%y/%m/%d},{}
)",
               fmt::format("{}", instance), fmt::format("{}", stat),
               tree->get_nb_nodes_explored(), fmt::localtime(result),
               fmt::format("{}", parms));
}

int Problem::to_screen() {
    int val = 0;

    switch (status) {
        case no_sol:
            fmt::print(
                "We didn't decide if this instance is feasible or "
                "infeasible\n");
            break;

        case feasible:
        case lp_feasible:
        case meta_heuristic:
            fmt::print(
                "A suboptimal schedule with relative error {} is found.\n",
                static_cast<double>(stat.global_upper_bound -
                                    stat.global_lower_bound) /
                    (stat.global_lower_bound));
            break;

        case optimal:
            fmt::print("The optimal schedule is found.\n");
            break;
    }

    fmt::print(
        "Compute_schedule took {} seconds(tot_scatter_search {}, "
        "tot_branch_and_bound {}, tot_lb_lp_root {}, tot_lb_lp {}, tot_lb "
        "{}, "
        "tot_pricing {}, tot_build_dd {}) and {} seconds in real time\n",
        stat.total_time_str(Statistics::cputime_timer, 2),
        stat.total_time_str(Statistics::heuristic_timer, 2),
        stat.total_time_str(Statistics::bb_timer, 2),
        stat.total_time_str(Statistics::lb_root_timer, 2),
        stat.total_time_str(Statistics::lb_timer, 2),
        stat.total_time_str(Statistics::lb_timer, 2),
        stat.total_time_str(Statistics::pricing_timer, 2),
        stat.total_time_str(Statistics::build_dd_timer, 2),
        stat.real_time_total);
    return val;
}
