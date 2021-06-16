#include "BranchNode.hpp"
#include <fmt/core.h>
#include <algorithm>
#include <memory>
#include <numeric>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/all.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/partial_sum.hpp>
#include <vector>
#include "Parms.h"
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "branch-and-bound/btree.h"
#include "util.h"
#include "wctprivate.h"

BranchNodeBase::BranchNodeBase(std::unique_ptr<NodeData> _pd, bool _isRoot)
    : State(_isRoot),
      pd(std::move(_pd)) {
    if (_isRoot) {
        set_ub(double(pd->opt_sol.tw));
        pd->build_rmp();
        pd->solve_relaxation();
        pd->stat.start_resume_timer(Statistics::lb_root_timer);
        pd->compute_lower_bound();
        pd->stat.suspend_timer(Statistics::lb_root_timer);
        set_lb(pd->lower_bound);
        if (pd->solver->get_is_integer_solution()) {
            set_obj_value(pd->LP_lower_bound);
        }
    }
}

void BranchNodeBase::branch(BTree* bt) {
    auto*       solver = get_pricersolver();
    const auto& instance = get_instance_info();
    const auto& parms = pd->parms;

    if (!parms.strong_branching && dbg_lvl() > 0) {
        fmt::print("\nDOING STRONG BRANCHING...\n\n");
    }

    if (bt->getGlobalUB() < pd->upper_bound) {
        pd->upper_bound = static_cast<int>(bt->getGlobalUB());
        solver->UB = bt->getGlobalUB();
    }

    auto& x_job_time = solver->calculate_job_time();

    // initialize the middle times and scores used for branching
    std::vector<BranchCand> best_cand{};

    for (auto&& [x_j, job] : ranges::views::zip(x_job_time, instance.jobs)) {
        auto aux_vec = ranges::views::iota(0UL, x_j.size()) |
                       ranges::views::transform([&](const auto& tmp) {
                           return job->weighted_tardiness_start(tmp);
                       }) |
                       ranges::to_vector;
        auto sum = ranges::inner_product(x_j, aux_vec, 0.0);
        auto lb_it = ranges::lower_bound(aux_vec, sum);
        if (*lb_it - sum > EPS &&
            std::min(sum - std::floor(sum), std::ceil(sum) - sum) > EPS) {
            auto tmp_t = ranges::distance(aux_vec.begin(), lb_it) - 1;
            best_cand.emplace_back(
                std::min(*lb_it - sum,
                         sum - job->weighted_tardiness_start(tmp_t)),
                job->job, tmp_t);
        }
    }

    if (best_cand.empty()) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto aux_vec =
                ranges::views::iota(0UL, x_j.size()) |
                ranges::views::transform([&](const auto& tmp) { return tmp; }) |
                ranges::to_vector;
            auto sum = ranges::inner_product(x_j, aux_vec, 0.0);
            auto lb_it = ranges::lower_bound(aux_vec, sum);
            if (*lb_it - sum > EPS &&
                std::min(sum - std::floor(sum), std::ceil(sum) - sum) > EPS) {
                auto tmp_t = ranges::distance(aux_vec.begin(), lb_it) - 1;
                best_cand.emplace_back(std::min(*lb_it - sum, sum - (tmp_t)),
                                       job->job, tmp_t);
            }
        }
    }

    if (best_cand.empty()) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto part_sum = ranges::views::partial_sum(x_j, std::plus<>{});
            auto br_point_it =
                ranges::lower_bound(part_sum, parms.branching_point - EPS);
            auto aux = std::min(*br_point_it - std::floor(*br_point_it),
                                std::ceil(*br_point_it) - *br_point_it);
            if (aux > EPS) {
                auto tmp_t = ranges::distance(part_sum.begin(), br_point_it);
                best_cand.emplace_back(aux, job->job, tmp_t);
            }
        }
    }

    if (best_cand.empty()) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto part_sum = ranges::views::partial_sum(x_j, std::plus<>{});
            auto br_point_it = ranges::upper_bound(part_sum, 0.0);
            auto aux = std::min(*br_point_it - std::floor(*br_point_it),
                                std::ceil(*br_point_it) - *br_point_it);
            if (aux > EPS) {
                auto tmp_t = ranges::distance(part_sum.begin(), br_point_it);
                best_cand.emplace_back(aux, job->job, tmp_t);
            }
        }
    }

    if (best_cand.empty()) {
        fmt::print(stderr, "ERROR: no branching found!\n");
        for (auto&& [j, job] : instance.jobs | ranges::views::enumerate) {
            fmt::print(stderr, "j={}:", j);
            for (auto&& [t, x] :
                 x_job_time[j] | ranges::views::enumerate |
                     ranges::views::filter(
                         [&](const auto& tmp) { return (tmp.second > EPS); })) {
                fmt::print(stderr, " ({},{})", t, x);
            }
            fmt::print(stderr, "\n");
        }

        exit(-1);
    }

    auto best_min_gain = 0.0;
    auto best_job = -1;
    auto best_time = 0;

    std::unique_ptr<BranchNodeBase> best_right = nullptr;
    std::unique_ptr<BranchNodeBase> best_left = nullptr;
    auto nb_cand = std::min(std::max(parms.strong_branching, 1),
                            static_cast<int>(best_cand.size()));
    best_cand |= ranges::actions::sort(std::greater<>{}, &BranchCand::score);

    auto nb_non_improvements = 0UL;
    for (auto& it : best_cand | ranges::views::take(nb_cand)) {
        auto data_nodes{pd->create_child_nodes(it.job, it.t)};
        std::array<std::unique_ptr<BranchNodeBase>, 2> child_nodes;
        std::array<double, 2>                          scores{};
        std::array<bool, 2>                            fathom{};
        auto                                           left = true;

        for (auto&& [data, node, score, f] :
             ranges::views::zip(data_nodes, child_nodes, scores, fathom)) {
            node = std::make_unique<BranchNodeBase>(std::move(data));
            node->compute_bounds(bt);

            auto cost = node->pd->LP_lower_bound;
            if (dbg_lvl() > 0) {
                fmt::print(
                    "STRONG BRANCHING {} PROBE: j = {}, t = {},"
                    " DWM LB = {:9.2f} in iterations {} \n\n",
                    left ? "LEFT" : "RIGHT", it.job, it.t,
                    cost + pd->instance.off, node->pd->iterations);
            }
            if (cost >= pd->opt_sol.tw - 1.0 + IntegerTolerance ||
                node->get_data_ptr()->solver->get_is_integer_solution()) {
                f = true;
            }
            score = (parms.scoring_parameter == weighted_sum_scoring_parameter)
                        ? cost / pd->LP_lower_bound
                        : std::abs(cost - pd->LP_lower_bound);
            left = false;
        }

        auto best_score = parms.scoring_function(scores[0], scores[1]);

        if (best_score > best_min_gain ||
            ranges::any_of(fathom, ranges::identity{})) {
            best_min_gain = best_score;
            best_right = std::move(child_nodes[1]);
            best_left = std::move(child_nodes[0]);
            nb_non_improvements = 0;
        } else {
            nb_non_improvements++;
        }

        if (ranges::any_of(fathom, std::identity{}) ||
            nb_non_improvements > 5) {
            break;
        }
    }

    /** Process the branching nodes insert them in the tree */
    bt->set_state_bounds_computed(true);
    bt->process_state(std::move(best_right));
    bt->process_state(std::move(best_left));
    bt->set_state_bounds_computed(false);
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
    set_feasible(pd->status != NodeData::infeasible);
    return (solver->get_is_integer_solution() || !is_feasible());
}

void BranchNodeBase::apply_final_pruning_tests(BTree* bt) {
    auto* pricing_solver = get_pricersolver();
    auto* pd_ptr = get_data_ptr();
    if (pricing_solver->get_nb_vertices() /
            double(pd->stat.size_graph_after_reduced_cost_fixing) <
        0.4) {
        fmt::print("MIPPING !!!!!\n");
        pricing_solver->build_mip();
    }
}

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

BranchCand::BranchCand(double _score, int _job, int _t)
    : score(_score),
      job(_job),
      t(_t) {}

BranchNodeRelBranching::BranchNodeRelBranching(std::unique_ptr<NodeData> _data,
                                               bool                      isRoot)
    : BranchNodeBase(std::move(_data), isRoot) {}

void BranchNodeRelBranching::branch(BTree* bt) {
    auto*       node_data = get_data_ptr();
    auto*       pricing_solver = get_pricersolver();
    const auto& instance = get_instance_info();

    if (bt->getGlobalUB() < node_data->upper_bound) {
        node_data->upper_bound = static_cast<int>(bt->getGlobalUB());
        pricing_solver->UB = bt->getGlobalUB();
    }

    auto& x_job_time = pricing_solver->calculate_job_time();

    std::vector<BranchCand> best_cand{};

    for (auto&& [x_j, job] : ranges::views::zip(x_job_time, instance.jobs)) {
        auto prev = -1;

        auto sum = ranges::inner_product(x_j, ranges::views::iota(0), 0.0);
        auto lb_sum_t = ranges::lower_bound(ranges::views::iota(0), sum);
        // for (auto&& [t, x] : x_j | ranges::views::enumerate) {
        //     if (t > sum) {
        //         break;
        //     }

        //     prev = t;
        // }

        if (true) {
            // best_cand.emplace_back(sum -, job->job, prev);
        }
    }

    std::unique_ptr<BranchNodeRelBranching> best_right = nullptr;
    std::unique_ptr<BranchNodeRelBranching> best_left = nullptr;

    auto best_min_gain = 0.0;
    auto best_job = -1;
    auto best_time = 0;

    bt->set_state_bounds_computed(true);
    bt->process_state(std::move(best_right));
    bt->process_state(std::move(best_left));
    bt->set_state_bounds_computed(false);
    bt->update_global_lb();

    if (dbg_lvl() > 1) {
        fmt::print("Branching choice: j = {}, t = {}, best_gain = {}\n",
                   best_job, best_time, best_min_gain + instance.off);
    }
}

const Instance& BranchNodeBase::get_instance_info() const {
    return pd->instance;
}

PricerSolverBase* BranchNodeBase::get_pricersolver() const {
    return pd->solver.get();
}
