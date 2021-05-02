#include "BranchNode.hpp"
#include <fmt/core.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <range/v3/action/sort.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/numeric/iota.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/take.hpp>
#include <span>
#include <vector>
#include "Parms.h"
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "branch-and-bound/btree.h"
#include "util.h"
#include "wctprivate.h"

using std::ranges::sort;

BranchNodeBase::BranchNodeBase(std::unique_ptr<NodeData> _pd, bool _isRoot)
    : State(_isRoot),
      pd(std::move(_pd)) {
    if (_isRoot) {
        pd->build_rmp();
        pd->solve_relaxation();
        pd->stat.start_resume_timer(Statistics::lb_root_timer);
        pd->compute_lower_bound();
        pd->stat.suspend_timer(Statistics::lb_root_timer);
        set_lb(pd->lower_bound);
        set_obj_value(pd->LP_lower_bound);
    }
}

void BranchNodeBase::branch(BTree* bt) {
    auto*       solver = pd->solver.get();
    auto        nb_jobs = pd->nb_jobs;
    auto        strong_branching = false;
    const auto& parms = pd->parms;

    if (!strong_branching && dbg_lvl() > 0) {
        fmt::print("\nDOING STRONG BRANCHING...\n\n");
    }

    auto fathom_left = false;
    auto fathom_right = false;
    if (bt->getGlobalUB() < pd->upper_bound) {
        pd->upper_bound = static_cast<int>(bt->getGlobalUB());
        solver->UB = bt->getGlobalUB();
    }

    // calculate the fraction of each job finishing at each time in the
    // relaxation solution
    std::vector<std::vector<double>> x_job_time(nb_jobs);
    for (auto& it : x_job_time) {
        it.resize(pd->instance.H_max + 1, 0.0);
    }

    solver->calculate_job_time(&x_job_time);

    // initialize the order for evaluating the branching jobs
    std::vector<size_t> ord(nb_jobs);
    ranges::iota(ord, 0UL);

    // initialize the middle times and scores used for branching
    std::vector<int>        middle_time(nb_jobs, -1);
    std::vector<BranchCand> best_cand(NumStrBrCandidates, BranchCand());
    std::vector<double>     branch_scores(nb_jobs, 0.0);
    std::vector<double>     min_completion_time(nb_jobs,
                                            std::numeric_limits<double>::max());
    std::vector<int>        lb_C(nb_jobs, std::numeric_limits<int>::max());
    std::vector<int>        ub_C(nb_jobs, 0);

    for (auto& i : ord) {
        auto  prev = -1;
        auto  accum = 0.0;
        auto  dist_zero = 0.0;
        auto* job = (solver->jobs)[i].get();
        for (auto&& [t, x] : x_job_time[i] | ranges::views::enumerate |
                                 ranges::views::filter([&](auto tmp) {
                                     return (tmp.second > EPS);
                                 })) {
            accum += x;
            if ((accum >= (1.0 - TargetBrTimeValue)) && (prev != -1) &&
                (middle_time[i] == -1)) {
                middle_time[i] = (t + job->processing_time + prev) / 2;
                branch_scores[i] =
                    double(job->weighted_tardiness(middle_time[i])) * accum -
                    dist_zero;
            }

            if (middle_time[i] != -1) {
                branch_scores[i] +=
                    double(job->weighted_tardiness(t + job->processing_time) -
                           job->weighted_tardiness(middle_time[i])) *
                    x;
            }

            dist_zero +=
                double(job->weighted_tardiness(t + job->processing_time)) * x;

            prev = t + job->processing_time;
            // accum += x;
            // if ((accum >= (1.0 - TargetBrTimeValue)) &&
            //     (middle_time[i] == -1)) {
            //     middle_time[i] = (t + job->processing_time);
            //     branch_scores[i] = accum;
            //     break;
            //     // double(job->weighted_tardiness(middle_time[i])) * accum -
            //     // dist_zero;
            // }

            // lb_C[i] =
            //     std::min(lb_C[i], static_cast<int>(t) +
            //     job->processing_time);

            // ub_C[i] =
            // std::max(ub_C[i], static_cast<int>(t) + job->processing_time);
            // if (middle_time[i] != -1) {
            //     branch_scores[i] +=
            //         double(job->weighted_tardiness(t + job->processing_time)
            //         -
            //                job->weighted_tardiness(middle_time[i])) *
            //         x;
            // }

            // dist_zero +=
            //     double(job->weighted_tardiness(t + job->processing_time)) *
            //     x;
        }

        // branch_scores[i] = avg_completion_time[i] - min_completion_time[i];
        // middle_time[i] = min_completion_time[i];

        // if (lb_C[i] != ub_C[i]) {
        //     fmt::print("test {}\n", (lb_C[i] + ub_C[i]) / 2);
        //     best_cand.emplace_back(i, (lb_C[i] + ub_C[i]) / 2, pd.get());
        // }

        auto minimum = std::min_element(best_cand.begin(), best_cand.end());

        if (minimum->score < branch_scores[i]) {
            minimum->score = branch_scores[i];
            minimum->job = i;
        }
    }

    auto best_min_gain = 0.0;
    auto best_job = -1;
    auto best_time = 0;

    std::unique_ptr<BranchNodeBase> best_right = nullptr;
    std::unique_ptr<BranchNodeBase> best_left = nullptr;
    best_cand |= ranges::actions::sort(std::less<BranchCand>());
    for (auto& it : best_cand | ranges::views::filter(
                                    [](auto& tmp) { return tmp.job != -1; })) {
        auto i = static_cast<size_t>(it.job);

        auto  left_gain = 0.0;
        auto  right_gain = 0.0;
        auto* job = pd->instance.jobs[i].get();
        for (auto&& [t, x] : x_job_time[i] | ranges::views::enumerate) {
            if (t + job->processing_time <= middle_time[i]) {
                left_gain += x;
            } else {
                right_gain += x;
            }
        }

        if (strong_branching) {
            auto min_gain = std::min(right_gain, left_gain);

            if (min_gain > best_min_gain) {
                best_min_gain = min_gain;
                best_job = i;
                best_time = middle_time[i];
            }
        } else {
            auto  left = pd->clone();
            auto* left_solver = left->solver.get();
            auto  left_node_branch =
                std::make_unique<BranchNodeBase>(std::move(left));
            left_solver->split_job_time(i, middle_time[i], true);
            left_node_branch->compute_bounds(bt);

            auto approx = left_gain;
            left_gain = left_node_branch->pd->LP_lower_bound;
            if (dbg_lvl() > 0) {
                fmt::print(
                    "STRONG BRANCHING LEFT PROBE: j = {}, t = {},"
                    " DWM LB = {:9.2f} in iterations {} ({})\n\n",
                    i, middle_time[i], left_gain + pd->instance.off,
                    left_node_branch->pd->iterations, approx);
            }
            if (left_gain >= pd->opt_sol.tw - 1.0 + IntegerTolerance ||
                left_node_branch->get_data_ptr()
                    ->solver->get_is_integer_solution()) {
                fathom_left = true;
            }
            left_gain =
                (parms.scoring_parameter == weighted_sum_scoring_parameter)
                    ? left_node_branch->pd->LP_lower_bound / pd->LP_lower_bound
                    : std::abs(left_node_branch->pd->LP_lower_bound -
                               pd->LP_lower_bound);

            // build the right node and solve its root LP only

            auto  right = pd->clone();
            auto* right_solver = right->solver.get();
            right_solver->split_job_time(i, middle_time[i], false);
            auto right_node_branch =
                std::make_unique<BranchNodeBase>(std::move(right));

            right_node_branch->compute_bounds(bt);

            approx = right_gain;
            right_gain = right_node_branch->pd->LP_lower_bound;
            if (dbg_lvl() > 0) {
                fmt::print(
                    "STRONG BRANCHING RIGHT PROBE: j = {}, t = {},"
                    " DWM LB = {:9.2f} in iterations {} ({})\n\n",
                    i, middle_time[i], right_gain + pd->instance.off,
                    right_node_branch->pd->iterations, approx);
            }
            if (right_gain >= pd->opt_sol.tw - 1.0 + IntegerTolerance ||
                right_node_branch->get_data_ptr()
                    ->solver->get_is_integer_solution()) {
                fathom_right = true;
            }

            right_gain =
                (parms.scoring_parameter == weighted_sum_scoring_parameter)
                    ? right_node_branch->pd->LP_lower_bound / pd->LP_lower_bound
                    : std::abs(right_node_branch->pd->LP_lower_bound -
                               pd->LP_lower_bound);

            // update the branching choice
            auto min_gain = parms.scoring_function(left_gain, right_gain);

            if (min_gain > best_min_gain) {
                best_min_gain = min_gain;
                best_job = i;
                best_time = middle_time[i];

                best_right = std::move(right_node_branch);
                best_left = std::move(left_node_branch);
            }

            if ((fathom_left || fathom_right)) {
                break;
            }
        }
    }

    if (best_cand.empty()) {
        fmt::print(stderr, "ERROR: no branching found!\n");
        for (auto& j : ord) {
            auto* job = (solver->jobs)[j].get();
            fmt::print(stderr, "j={}:", j);
            for (auto&& [t, x] : x_job_time[j] | ranges::views::enumerate) {
                if (x > ERROR) {
                    fmt::print(stderr, " ({},{})", t + job->processing_time, x);
                }
            }
            fmt::print(stderr, "\n");
        }

        exit(-1);
    }

    /** Process the branching nodes insert them in the tree */
    if (strong_branching) {
        auto  left = pd->clone();
        auto* left_solver = left->solver.get();
        auto  left_node_branch =
            std::make_unique<BranchNodeBase>(std::move(left));
        left_solver->split_job_time(best_job, best_time, false);
        left_node_branch->pd->branch_job = best_job;
        left_node_branch->pd->completiontime = best_time;
        left_node_branch->pd->less = 0;
        bt->process_state(std::move(left_node_branch));

        auto  right = pd->clone();
        auto* right_solver = right->solver.get();
        auto  right_node_branch =
            std::make_unique<BranchNodeBase>(std::move(right));
        right_solver->split_job_time(best_job, best_time, true);
        right_node_branch->pd->branch_job = best_job;
        right_node_branch->pd->completiontime = best_time;
        right_node_branch->pd->less = 1;
        bt->process_state(std::move(right_node_branch));
    } else {
        bt->set_state_bounds_computed(true);
        best_left->pd->branch_job = best_right->pd->branch_job = best_job;
        best_left->pd->completiontime = best_right->pd->completiontime =
            best_time;
        best_left->pd->less = 0;
        best_right->pd->less = 1;
        bt->process_state(std::move(best_right));
        bt->process_state(std::move(best_left));
        bt->set_state_bounds_computed(false);
    }

    bt->update_global_lb();

    if (dbg_lvl() > 1) {
        fmt::print("Branching choice: j = {}, t = {}, best_gain = {}\n",
                   best_job, best_time, best_min_gain + pd->instance.off);
    }
}

void BranchNodeBase::compute_bounds(BTree* bt) {
    pd->build_rmp();
    pd->solve_relaxation();
    pd->compute_lower_bound();
    set_lb(pd->lower_bound);
    if (pd->solver->get_is_integer_solution()) {
        set_obj_value(pd->lower_bound);
    }
}

void BranchNodeBase::assess_dominance([[maybe_unused]] State* otherState) {}

bool BranchNodeBase::is_terminal_state() {
    auto* solver = pd->solver.get();
    return solver->get_is_integer_solution();
}

void BranchNodeBase::apply_final_pruning_tests(BTree* bt) {}

void BranchNodeBase::update_data(double upper_bound) {
    // auto* statistics = pd->stat;
    pd->stat.global_upper_bound = static_cast<int>(upper_bound);
}

void BranchNodeBase::print(const BTree* bt) const {
    auto* solver = pd->solver.get();
    if (get_is_root_node()) {
        fmt::print("{0:^10}|{1:^30}|{2:^30}|{3:^10}|{4:^10}|\n", "Nodes",
                   "Current Node", "Objective Bounds", "Branch", "Work");
        fmt::print(
            R"({0:^5}{1:^5}|{2:^10}{3:^10}{10:^10}|{4:>10}{5:>10}{6:>10}|{7:>5}{8:>5}|{9:>5}{11:>5}
)",
            "Expl", "UnEx", "Obj", "Depth", "Primal", "Dual", "Gap", "Job",
            "Time", "Time", "Size", "Iter");
    }
    fmt::print(
        R"({0:>5}{1:>5}|{2:10.2f}{3:>10}{4:>10}|{5:10.2f}{6:10.2f}{7:10.2f}|{8:>5}{9:>5}|{10:>5}{11:>5}|
)",
        bt->get_nb_nodes_explored(), bt->get_nb_nodes(),
        pd->LP_lower_bound + pd->instance.off, pd->depth,
        solver->get_nb_vertices(), bt->getGlobalUB() + pd->instance.off,
        bt->getGlobalLB() + pd->instance.off, 0.0, pd->branch_job,
        pd->completiontime, bt->get_run_time_start(), pd->iterations);
}

BranchCand::BranchCand(int _job, int _t, const NodeData* parrent)
    : job(_job),
      t(_t) {
    //   left(std::move(std::make_unique<BranchNodeBase>(parrent->clone()))),
    //   right(std::move(std::make_unique<BranchNodeBase>(parrent->clone()))) {
    // auto* left_data_ptr = left->get_data_ptr();
    // auto* right_data_ptr = right->get_data_ptr();

    // left_data_ptr->solver->split_job_time(job, t, true);
    // left_data_ptr->build_rmp();
    // left_data_ptr->solve_relaxation();
    // left_data_ptr->estimate_lower_bound(10);

    // right_data_ptr->solver->split_job_time(job, t, false);
    // right_data_ptr->build_rmp();
    // right_data_ptr->solve_relaxation();
    // right_data_ptr->estimate_lower_bound(10);

    // score = std::min(right_data_ptr->LP_lower_bound -
    // parrent->LP_lower_bound,
    //                  left_data_ptr->LP_lower_bound -
    //                  parrent->LP_lower_bound);
}
