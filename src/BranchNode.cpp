#include "BranchNode.hpp"
#include <fmt/core.h>
#include <algorithm>
#include <limits>
#include <memory>
#include <span>
#include "PricerSolverBase.hpp"
#include "branch-and-bound/btree.h"
#include "solver.h"
#include "util.h"
#include "wctprivate.h"

BranchNodeBase::BranchNodeBase(std::unique_ptr<NodeData> _pd, bool _isRoot)
    : State(_isRoot),
      pd(std::move(_pd)) {
    if (_isRoot) {
        pd->build_rmp();
        pd->solve_relaxation();
        pd->compute_lower_bound();
        set_lb(pd->lower_bound);
        set_obj_value(pd->LP_lower_bound);
    }
}

void BranchNodeBase::branch(BTree* bt) {
    auto* solver = pd->solver.get();
    auto  nb_jobs = pd->nb_jobs;
    auto  strong_branching =
        ((pd->opt_sol.tw - pd->LP_lower_bound) < (1.0 + IntegerTolerance));

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
    for (auto i = 0; i < nb_jobs; i++) {
        x_job_time[i].resize(pd->instance.H_max + 1, 0.0);
    }

    solver->calculate_job_time(&x_job_time);

    // initialize the order for evaluating the branching jobs
    std::vector<int> ord(nb_jobs);
    for (int k = 0; k < nb_jobs; k++) {
        ord[k] = k;
    }

    // initialize the middle times and scores used for branching
    std::vector<int>        middle_time(nb_jobs, -1);
    std::vector<BranchCand> best_cand(NumStrBrCandidates, BranchCand());
    std::vector<double>     branch_scores(nb_jobs, 0.0);
    std::vector<double>     avg_completion_time(nb_jobs, 0.0);
    std::vector<double>     min_completion_time(nb_jobs,
                                            std::numeric_limits<double>::max());

    for (auto k = 0; k < nb_jobs; k++) {
        auto  i = ord[k];
        auto  prev = -1;
        auto  accum = 0.0;
        auto  dist_zero = 0.0;
        auto* job = (solver->jobs)[i].get();
        for (auto t = 0; t < pd->instance.H_max + 1; t++) {
            accum += x_job_time[i][t];
            // avg_completion_time[i] +=
            //     x_job_time[i][t] * (t + job->processing_time);
            // if (x_job_time[i][t] > 1e-5) {
            //     min_completion_time[i] = std::min(
            //         min_completion_time[i], double(t +
            //         job->processing_time));
            // }

            if ((accum >= (1.0 - TargetBrTimeValue)) &&
                (x_job_time[i][t] > EPS) && (prev != -1) &&
                (middle_time[i] == -1)) {
                middle_time[i] = (t + job->processing_time + prev) / 2;
                branch_scores[i] = double(middle_time[i]) * accum - dist_zero;
            }

            if (middle_time[i] != -1) {
                branch_scores[i] +=
                    double(t + job->processing_time - middle_time[i] + 1) *
                    x_job_time[i][t];
            }

            dist_zero += double(t + job->processing_time) * x_job_time[i][t];

            if (x_job_time[i][t] > EPS) {
                prev = t + job->processing_time;
            }
        }

        // branch_scores[i] = avg_completion_time[i] - min_completion_time[i];
        // middle_time[i] = min_completion_time[i];

        auto minimum = std::min_element(best_cand.begin(), best_cand.end());

        if (minimum->score < branch_scores[i]) {
            minimum->score = branch_scores[i];
            minimum->job = i;
        }
    }

    auto best_min_gain = 0.0;
    auto best_job = -1;
    auto best_time = 0;
    // std::span jobsarray{pd->jobarray->pdata, pd->jobarray->len};
    std::unique_ptr<BranchNodeBase> best_right = nullptr;
    std::unique_ptr<BranchNodeBase> best_left = nullptr;
    for (auto& it : best_cand) {
        auto i = it.job;
        if (i == -1) {
            continue;
        }

        auto  left_gain = 0.0;
        auto  right_gain = 0.0;
        auto* job = pd->instance.jobs[i].get();
        for (auto t = 0; t < pd->instance.H_max + 1; t++) {
            if (t + job->processing_time <= middle_time[i]) {
                left_gain += x_job_time[i][t];
            } else {
                right_gain += x_job_time[i][t];
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
            // std::unique_ptr<BranchNodeBase>(new BranchNodeBase(left));
            left_solver->split_job_time(i, middle_time[i], true);
            left_node_branch->compute_bounds(bt);

            auto approx = left_gain;
            left_gain = left_node_branch->get_obj_value();
            if (dbg_lvl() > 0) {
                fmt::print(
                    "STRONG BRANCHING LEFT PROBE: j = {}, t = {},"
                    " DWM LB = {:9.2f} in iterations {} ({})\n\n",
                    i, middle_time[i], left_gain + pd->instance.off,
                    left_node_branch->pd->iterations, approx);
            }
            if (left_gain >= pd->opt_sol.tw - 1.0 + IntegerTolerance ||
                left_solver->get_is_integer_solution()) {
                fathom_left = true;
            }

            // build the right node and solve its root LP only

            auto  right = pd->clone();
            auto* right_solver = right->solver.get();
            auto  right_node_branch =
                std::make_unique<BranchNodeBase>(std::move(right));

            right_solver->split_job_time(i, middle_time[i], false);
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
                right_solver->get_is_integer_solution()) {
                fathom_right = true;
            }

            // update the branching choice
            auto min_gain = std::min(right_gain, left_gain);

            if (min_gain > best_min_gain || fathom_left || fathom_right) {
                best_min_gain = min_gain;
                best_job = i;
                best_time = middle_time[i];

                best_right = std::move(right_node_branch);
                best_left = std::move(left_node_branch);
            }

            if (fathom_left || fathom_right) {
                break;
            }
        }
    }

    if (best_job == -1) {
        fmt::print(stderr, "ERROR: no branching found!\n");
        for (auto k = 0; k < nb_jobs; k++) {
            auto  j = ord[k];
            auto* job = (solver->jobs)[j].get();
            fmt::print(stderr, "j={}:", j);
            for (int t = 0; t < pd->instance.H_max; t++) {
                if (x_job_time[j][t] > ERROR) {
                    fmt::print(stderr, " ({},{})", t + job->processing_time,
                               x_job_time[j][t]);
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
        bt->process_state(std::move(best_left));
        bt->process_state(std::move(best_right));
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
    set_obj_value(pd->LP_lower_bound);
}

void BranchNodeBase::assess_dominance(State* otherState) {}

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
