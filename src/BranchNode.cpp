#include "BranchNode.hpp"
#include <fmt/core.h>
#include <fmt/format.h>
#include <algorithm>
#include <array>
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <range/v3/action/shuffle.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/heap_algorithm.hpp>
#include <range/v3/algorithm/lower_bound.hpp>
#include <range/v3/all.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/partial_sum.hpp>
#include <range/v3/view/zip.hpp>
#include <vector>
#include "NodeData.h"
#include "Parms.h"
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "branch-and-bound/btree.h"
#include "orutils/util.h"

namespace fmt {
template <>
struct formatter<boost::multiprecision::cpp_int> {
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx) {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format([[maybe_unused]] const boost::multiprecision::cpp_int& p,
                FormatContext&                                         ctx) {
        return format_to(ctx.begin(), "boost::number");
    }
};
}  // namespace fmt

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
    // auto&       stat = pd->stat;

    if (!parms.strong_branching && dbg_lvl() > 0) {
        fmt::print("\nDOING STRONG BRANCHING...\n\n");
    }

    if (bt->getGlobalUB() < pd->upper_bound) {
        pd->upper_bound = static_cast<int>(bt->getGlobalUB());
        solver->UB = bt->getGlobalUB();
    }

    auto& x_job_time = solver->calculate_job_time();
    auto  x_ref = solver->get_pair_x();

    std::array<BranchCandidate, NbCandidates> candidates{};

    ranges::make_heap(candidates, std::greater<>{},
                      [](const auto& tmp) { return tmp.score; });

    auto tmp_t = 0;
    auto rng_zip = ranges::views::zip(x_job_time, instance.jobs);
    auto rng_sum_wt =
        rng_zip | ranges::views::filter([&tmp_t](const auto& tmp) {
            auto aux_vec = ranges::views::iota(size_t{}, tmp.first.size()) |
                           ranges::views::transform([&](const auto& t) {
                               return tmp.second->weighted_tardiness_start(t);
                           }) |
                           ranges::to<std::vector>;
            auto sum = ranges::inner_product(tmp.first, aux_vec, 0.0);
            auto lb_it = ranges::lower_bound(aux_vec, sum);
            tmp_t = ranges::distance(aux_vec.begin(), lb_it) - 1;
            return std::min(sum - std::floor(sum), std::ceil(sum) - sum) > EPS;
        });

    for (auto&& [x_j, job] : rng_sum_wt | ranges::views::take(std::min(
                                              size_t{50}, instance.nb_jobs))) {
        while ((x_ref.second[job->job][tmp_t] ==
                x_ref.second[job->job][tmp_t + 1]) &&
               job->weighted_tardiness_start(tmp_t) == 0) {
            tmp_t += 1;
        }
        auto aux = BranchCandidate(pd->get_score_value(),
                                   pd->create_child_nodes(job->job, tmp_t));
        ranges::make_heap(candidates, std::greater<>{},
                          [](const auto& tmp) { return tmp.score; });
        if (aux.score > candidates.front().score) {
            std::swap(aux, candidates.front());
        }
        // }
    }

    if (ranges::all_of(candidates, [](const auto& tmp) { return tmp.empty; })) {
        auto rng_t =
            rng_zip | ranges::views::filter([&tmp_t](const auto& tmp) {
                auto aux_vec = ranges::views::iota(0UL, tmp.first.size()) |
                               ranges::to<std::vector>;
                auto sum = ranges::inner_product(tmp.first, aux_vec, 0.0);
                auto lb_it = ranges::lower_bound(aux_vec, sum);
                tmp_t = ranges::distance(aux_vec.begin(), lb_it) - 1;
                return std::min(sum - std::floor(sum), std::ceil(sum) - sum) >
                       EPS;
            });
        for (auto&& [x_j, job] : rng_t | ranges::views::take(std::min(
                                             size_t{50}, instance.nb_jobs))) {
            while ((x_ref.second[job->job][tmp_t] ==
                    x_ref.second[job->job][tmp_t + 1]) &&
                   (job->weighted_tardiness_start(tmp_t) == 0)) {
                tmp_t += 1;
            }
            auto aux = BranchCandidate(pd->get_score_value(),
                                       pd->create_child_nodes(job->job, tmp_t));
            ranges::make_heap(candidates, std::greater<>{},
                              [](const auto& tmp) { return tmp.score; });
            if (aux.score > candidates.front().score) {
                std::swap(aux, candidates.front());
            }
            // }
        }
    }

    if (ranges::all_of(candidates, [](const auto& tmp) { return tmp.empty; })) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto part_sum = ranges::views::partial_sum(x_j, std::plus<>{});
            auto br_point_it =
                ranges::lower_bound(part_sum, parms.branching_point - EPS);
            auto aux = std::min(*br_point_it - std::floor(*br_point_it),
                                std::ceil(*br_point_it) - *br_point_it);
            if (aux > EPS) {
                tmp_t = ranges::distance(part_sum.begin(), br_point_it);
                while ((x_ref.second[job->job][tmp_t] ==
                        x_ref.second[job->job][tmp_t + 1]) &&
                       (job->weighted_tardiness_start(tmp_t))) {
                    tmp_t += 1;
                }
                auto aux_data =
                    BranchCandidate(pd->get_score_value(),
                                    pd->create_child_nodes(job->job, tmp_t));
                ranges::make_heap(candidates, std::greater<>{},
                                  [](const auto& tmp) { return tmp.score; });
                if (aux_data.score > candidates.front().score) {
                    std::swap(aux_data, candidates.front());
                }
            }
        }
    }

    if (ranges::all_of(candidates, [](const auto& tmp) { return tmp.empty; })) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto part_sum = ranges::views::partial_sum(x_j, std::plus<>{});
            auto br_point_it = ranges::upper_bound(part_sum, 0.0);
            auto aux = std::min(*br_point_it - std::floor(*br_point_it),
                                std::ceil(*br_point_it) - *br_point_it);
            if (aux > EPS) {
                tmp_t = ranges::distance(part_sum.begin(), br_point_it);
                while ((x_ref.second[job->job][tmp_t] ==
                        x_ref.second[job->job][tmp_t + 1]) &&
                       (job->weighted_tardiness_start(tmp_t))) {
                    tmp_t += 1;
                }
                auto aux_data =
                    BranchCandidate(pd->get_score_value(),
                                    pd->create_child_nodes(job->job, tmp_t));
                ranges::make_heap(candidates, std::greater<>{},
                                  [](const auto& tmp) { return tmp.score; });
                if (aux_data.score > candidates.front().score) {
                    std::swap(aux_data, candidates.front());
                }
            }
        }
    }

    if (ranges::all_of(candidates, [](const auto& tmp) { return tmp.empty; })) {
        fmt::print(stderr, "ERROR: no branching found!\n");
        for (auto&& [j, job] : instance.jobs | ranges::views::enumerate) {
            fmt::print(stderr, "j={}:", j);
            auto rng = x_job_time[j] | ranges::views::enumerate |
                       ranges::views::filter(
                           [&](const auto& tmp) { return (tmp.second > EPS); });
            for (auto&& [t, x] : rng) {
                fmt::print(stderr, " ({},{})", t, x);
            }
            fmt::print(stderr, "\n");
        }
        exit(-1);
    }

    auto best_min_gain = std::numeric_limits<double>::min();
    auto best_job = -1;
    auto best_time = 0;

    std::array<std::unique_ptr<BranchNodeBase>, 2> best{};
    auto nb_cand = std::min(std::max(parms.strong_branching, 1),
                            static_cast<int>(candidates.size()));

    candidates |= ranges::actions::sort(
        std::greater<>{}, [](const auto& tmp) { return tmp.score; });

    auto rng_best_cand =
        candidates |
        ranges::views::filter([](const auto& tmp) { return !(tmp.empty); });

    auto nb_non_improvements = 0UL;
    for (auto& it : rng_best_cand | ranges::views::take(nb_cand)) {
        std::array<std::unique_ptr<BranchNodeBase>, 2> child_nodes;
        std::array<double, 2>                          scores{};
        std::array<bool, 2>                            fathom{};
        auto                                           left = false;

        for (auto&& [data, node, score, f] : ranges::views::zip(
                 it.data_child_nodes, child_nodes, scores, fathom)) {
            node = std::make_unique<BranchNodeBase>(std::move(data));
            node->compute_bounds(bt);

            auto cost = node->pd->LP_lower_bound;
            if (dbg_lvl() > 0) {
                fmt::print(
                    "STRONG BRANCHING {} PROBE: j = {}, t = {},"
                    " DWM LB = {:9.2f} in iterations {} \n\n",
                    left ? "LEFT" : "RIGHT", node->pd->branch_job,
                    node->pd->completiontime, cost + pd->instance.off,
                    node->pd->iterations);
            }
            if (cost >= bt->getGlobalUB() - 1.0 + IntegerTolerance ||
                node->get_data_ptr()->solver->get_is_integer_solution()) {
                f = true;
            }
            score = (parms.scoring_parameter == weighted_sum_scoring_parameter)
                        ? cost / pd->LP_lower_bound
                        : std::abs(cost - pd->LP_lower_bound);
            left = true;
        }

        auto best_score = parms.scoring_function(scores);

        if (best_score > best_min_gain ||
            ranges::any_of(fathom, ranges::identity{})) {
            best_min_gain = best_score;
            best = std::move(child_nodes);
            nb_non_improvements = 0;
        } else {
            ++nb_non_improvements;
        }

        if (ranges::any_of(fathom, std::identity{}) ||
            nb_non_improvements > 3) {
            break;
        }
    }

    /** Process the branching nodes insert them in the tree */
    bt->set_state_bounds_computed(true);
    for (auto& it : best) {
        bt->process_state(std::move(it));
    }
    bt->set_state_bounds_computed(false);
    bt->update_global_lb();

    if (dbg_lvl() > 1) {
        fmt::print("Branching choice: j = {}, t = {}, best_gain = {}\n",
                   best_job, best_time, best_min_gain + pd->instance.off);
    }
}

void BranchNodeBase::compute_bounds([[maybe_unused]] BTree* bt) {
    if (!pd->solver->structure_feasible()) {
        pd->status = NodeData::infeasible;
        pd->LP_lower_bound_dual = pd->LP_lower_bound = pd->LP_lower_bound_BB =
            pd->upper_bound;
        pd->lower_bound = pd->upper_bound;
        return;
    }
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

void BranchNodeBase::apply_final_pruning_tests([[maybe_unused]] BTree* bt) {
    auto* pricing_solver = get_pricersolver();
    // auto* pd_ptr = get_data_ptr();
    if (pricing_solver->print_num_paths() < MAX_NB_PATHS) {
        pricing_solver->project_sol_on_original_variables(pd->opt_sol);
        pricing_solver->build_mip();
        auto feasible_aux = pricing_solver->evaluate_mip_model();
        set_dominated(!feasible_aux);
    }
}

void BranchNodeBase::update_data(double upper_bound) {
    // auto* statistics = pd->stat;
    pd->stat.global_upper_bound = static_cast<int>(upper_bound);
}

void BranchNodeBase::print(const BTree* bt) const {
    auto* solver = pd->solver.get();
    if (get_is_root_node()) {
        fmt::print("{0:^10}|{1:^75}|{2:^30}|{3:^10}|{4:^10}|\n", "Nodes",
                   "Current Node", "Objective Bounds", "Branch", "Work");
        fmt::print(
            R"({0:^5}{1:^5}|{2:>10}{3:>10}{10:>10}{12:>45}|{4:>10}{5:>10}{6:>10}|{7:>5}{8:>5}|{9:>5}{11:>5}|
)",
            "Expl", "UnEx", "Obj", "Depth", "Primal", "Dual", "Gap", "Job",
            "Time", "Time", "Size", "Iter", "NB Paths");
    }
    fmt::print(
        R"({0:>5}{1:>5}|{2:10.2f}{3:>10}{4:>10}{12:>45}|{5:10.2f}{6:10.2f}{7:10.2f}|{8:>5}{9:>5}|{10:>5}{11:>5}|
)",
        bt->get_nb_nodes_explored(), bt->get_nb_nodes(),
        pd->LP_lower_bound + pd->instance.off, pd->depth,
        solver->get_nb_vertices(), bt->getGlobalUB() + pd->instance.off,
        bt->getGlobalLB() + pd->instance.off, 0.0, pd->branch_job,
        pd->completiontime, bt->get_run_time_start(), pd->iterations,
        solver->print_num_paths().str());
}

BranchCandidate::BranchCandidate(
    double                                     _score,
    std::array<std::unique_ptr<NodeData>, 2>&& child_nodes)
    : score(_score),
      empty(false),
      data_child_nodes(std::move(child_nodes)) {
    const auto&           parms = data_child_nodes[0]->parms;
    std::array<double, 2> scores{};

    for (auto&& [node, s] : ranges::views::zip(data_child_nodes, scores)) {
        s = (parms.scoring_parameter == weighted_sum_scoring_parameter)
                ? _score / node->get_score_value()
                : std::abs(_score - node->get_score_value());
    }
    score = parms.scoring_function(scores);
}

BranchCandidate::BranchCandidate(double _score) : score(_score) {}

BranchNodeRelBranching::BranchNodeRelBranching(std::unique_ptr<NodeData> _data,
                                               bool                      isRoot)
    : BranchNodeBase(std::move(_data), isRoot) {}

void BranchNodeRelBranching::branch(BTree* bt) {
    auto*       node_data = get_data_ptr();
    auto*       pricing_solver = get_pricersolver();
    auto&       parms = node_data->parms;
    const auto& instance = get_instance_info();
    const auto  Tmax = instance.H_max;

    if (bt->getGlobalUB() < node_data->upper_bound) {
        node_data->upper_bound = static_cast<int>(bt->getGlobalUB());
        pricing_solver->UB = bt->getGlobalUB();
    }

    std::array<BranchCandidate, NbCandidates> candidates{};
    auto&       x_job_time = pricing_solver->calculate_job_time();
    const auto& x_ref = pricing_solver->get_pair_x();
    auto&       aux_branch_history = bt->get_branch_history();

    for (auto& it : aux_branch_history) {
        auto j = it.first / (Tmax + 1);
        auto t = it.first % (Tmax + 1);
        auto aux_x = x_ref.second[j][t];
        if (std::min(aux_x - std::floor(aux_x), std::ceil(aux_x) - aux_x) >
            EPS) {
            auto aux = BranchCandidate(it.second.get_score());
            ranges::make_heap(candidates, std::greater<>{},
                              [](const auto& tmp) { return tmp.score; });
            if (aux.score > candidates.front().score) {
                std::swap(aux, candidates.front());
            }
        }
    }

    for (auto&& [z_j, job] : ranges::views::zip(x_ref.second, instance.jobs)) {
        for (auto&& [t, x_j_t] : z_j | ranges::views::enumerate) {
            auto left_gain = x_j_t - std::floor(x_j_t);
            auto right_gain = std::ceil(x_j_t) - x_j_t;
            if (std::min(left_gain, right_gain) > EPS) {
                auto aux_it =
                    aux_branch_history.find(job->job * instance.H_max + t);
                if (aux_it != aux_branch_history.end()) {
                    // candidates.emplace_back(
                    //     aux_it->second.compute_score(parms.scoring_function,
                    //                                  left_gain, right_gain),
                    //     node_data->create_child_nodes(job->job, t));
                }
            }
        }
    }

    for (auto&& [x_j, job] : ranges::views::zip(x_job_time, instance.jobs)) {
        auto aux_vec = ranges::views::iota(size_t{}, x_j.size()) |
                       ranges::views::transform([&](const auto& tmp) {
                           return job->weighted_tardiness_start(tmp);
                       }) |
                       ranges::to_vector;
        auto sum = ranges::inner_product(x_j, aux_vec, 0.0);
        auto lb_it = ranges::lower_bound(aux_vec, sum);
        if (*lb_it - sum > EPS &&
            std::min(sum - std::floor(sum), std::ceil(sum) - sum) > EPS) {
            // auto tmp_t = ranges::distance(aux_vec.begin(), lb_it) - 1;
            // candidates.emplace_back(
            //     node_data->LP_lower_bound,
            //     node_data->create_child_nodes(job->job, tmp_t));
        }
    }

    if (candidates.empty()) {
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
                // auto tmp_t = ranges::distance(aux_vec.begin(), lb_it) - 1;
                // candidates.emplace_back(
                //     node_data->LP_lower_bound,
                //     node_data->create_child_nodes(job->job, tmp_t));
            }
        }
    }

    if (candidates.empty()) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto part_sum = ranges::views::partial_sum(x_j, std::plus<>{});
            auto br_point_it =
                ranges::lower_bound(part_sum, parms.branching_point - EPS);
            auto aux = std::min(*br_point_it - std::floor(*br_point_it),
                                std::ceil(*br_point_it) - *br_point_it);
            if (aux > EPS) {
                // auto tmp_t = ranges::distance(part_sum.begin(), br_point_it);
                // candidates.emplace_back(
                //     node_data->LP_lower_bound,
                //     node_data->create_child_nodes(job->job, tmp_t));
            }
        }
    }

    if (candidates.empty()) {
        for (auto&& [x_j, job] :
             ranges::views::zip(x_job_time, instance.jobs)) {
            auto part_sum = ranges::views::partial_sum(x_j, std::plus<>{});
            auto br_point_it = ranges::upper_bound(part_sum, 0.0);
            auto aux = std::min(*br_point_it - std::floor(*br_point_it),
                                std::ceil(*br_point_it) - *br_point_it);
            if (aux > EPS) {
                // auto tmp_t = ranges::distance(part_sum.begin(), br_point_it);
                // candidates.emplace_back(
                //     node_data->LP_lower_bound,
                //     node_data->create_child_nodes(job->job, tmp_t));
            }
        }
    }

    candidates |=
        ranges::actions::sort(std::greater<>{}, &BranchCandidate::score);

    // for (auto& it : candidates) {
    // }

    std::array<std::unique_ptr<BranchNodeRelBranching>, 2> best = {nullptr,
                                                                   nullptr};
    // std::unique_ptr<BranchNodeRelBranching>                best_left =
    // nullptr;

    auto best_min_gain = 0.0;
    auto best_job = -1;
    auto best_time = 0;

    bt->set_state_bounds_computed(true);
    for (auto& it : best) {
        bt->process_state(std::move(it));
    }
    // bt->process_state(std::move(best_right));
    // bt->process_state(std::move(best_left));
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
