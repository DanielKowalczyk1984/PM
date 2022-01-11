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

#include "PricerSolverZddForward.hpp"
#include <fmt/core.h>                            // for print
#include <array>                                 // for array
#include <iostream>                              // for operator<<, basic_os...
#include <memory>                                // for __shared_ptr_access
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=, basic_it...
#include <range/v3/view/drop.hpp>                // for drop, drop_fn
#include <range/v3/view/join.hpp>                // for join_view, join_view...
#include <range/v3/view/subrange.hpp>            // for subrange
#include <range/v3/view/take.hpp>                // for take_view, take, tak...
#include <range/v3/view/view.hpp>                // for operator|, view_closure
#include <span>                                  // for span
#include <vector>                                // for vector
#include "Instance.h"                            // for Instance
#include "Job.h"                                 // for Job
#include "Label.hpp"                             // for Label
#include "ModernDD/NodeBddStructure.hpp"         // for DdStructure
#include "ModernDD/NodeBddTable.hpp"             // for NodeTableEntity, Tab...
#include "PricerSolverBase.hpp"                  // for PricerSolverBase::RC...
#include "PricerSolverZdd.hpp"                   // for PricerSolverZdd
#include "ZddNode.hpp"                           // for SubNodeZdd, NodeZdd
/**
 *  zdd solver pricersolver for the flow formulation
 */

PricerSolverSimple::PricerSolverSimple(const Instance& instance)
    : PricerSolverZdd(instance) {
    fmt::print("Constructing ZDD with Forward Simple evaluator\n");
    fmt::print("number vertices ZDD = {}\n", get_nb_vertices());
    fmt::print("number edges ZDD = {}\n", get_nb_edges());
    evaluator = ForwardZddSimpleDouble(convex_constr_id);
    reversed_evaluator = BackwardZddSimpleDouble(convex_constr_id);
}

auto PricerSolverSimple::pricing_algorithm(double* _pi) -> PricingSolution {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

auto PricerSolverSimple::pricing_algorithm(std::span<const double>& _pi)
    -> PricingSolution {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

void PricerSolverSimple::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

void PricerSolverSimple::compute_labels(std::span<const double>& _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

auto PricerSolverSimple::evaluate_nodes(std::span<const double>& pi) -> bool {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();

    nb_removed_edges = 0;

    // /** check for each node the Lagrangian dual */
    for (auto& it : table |
                        ranges::views::take(decision_diagram->topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::join) {
        for (auto& iter : it.list) {
            int    w = iter->get_weight();
            Job*   job = it.get_job();
            double result = iter->forward_label[0].get_f() +
                            iter->y->backward_label[0].get_f() -
                            job->weighted_tardiness_start(w) + pi[job->job];
            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result >
                    UB - 1 + RC_FIXING &&
                (iter->calc_yes)) {
                iter->calc_yes = false;
                nb_removed_edges++;
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);

    return nb_removed_edges;
}

auto PricerSolverSimple::evaluate_nodes(double* pi) -> bool {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();

    nb_removed_edges = 0;

    // /** check for each node the Lagrangian dual */
    for (auto& it : table |
                        ranges::views::take(decision_diagram->topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::join) {
        for (auto& iter : it.list) {
            int    w = iter->get_weight();
            Job*   job = it.get_job();
            double result = iter->forward_label[0].get_f() +
                            iter->y->backward_label[0].get_f() -
                            job->weighted_tardiness_start(w) + pi[job->job];
            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result >
                    UB - 1 + RC_FIXING &&
                (iter->calc_yes)) {
                iter->calc_yes = false;
                nb_removed_edges++;
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);

    return nb_removed_edges;
}

PricerSolverZddCycle::PricerSolverZddCycle(const Instance& instance)
    : PricerSolverZdd(instance) {
    fmt::print("Constructing ZDD with Forward ZddCycle evaluator\n");
    fmt::print("number vertices ZDD = {}", get_nb_vertices());
    fmt::print("number edges ZDD = {}", get_nb_edges());
    evaluator = ForwardZddCycleDouble(convex_constr_id);
    reversed_evaluator = BackwardZddCycleDouble(convex_constr_id);
}

auto PricerSolverZddCycle::pricing_algorithm(double* _pi) -> PricingSolution {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

auto PricerSolverZddCycle::pricing_algorithm(std::span<const double>& _pi)
    -> PricingSolution {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

void PricerSolverZddCycle::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

void PricerSolverZddCycle::compute_labels(std::span<const double>& _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

auto PricerSolverZddCycle::evaluate_nodes(std::span<const double>& pi) -> bool {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    auto reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (auto& it : table |
                        ranges::views::take(decision_diagram->topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::join) {
        auto* job = it.get_job();
        for (auto& iter : it.list) {
            auto w = iter->get_weight();

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (iter->forward_label[0].prev_job_forward() != job) {
                if (iter->y->backward_label[0].prev_job_backward() != job) {
                    auto result = iter->forward_label[0].get_f() +
                                  iter->y->backward_label[0].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    auto result = iter->forward_label[0].get_f() +
                                  iter->y->backward_label[1].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            } else {
                if (iter->y->backward_label[0].prev_job_backward() != job) {
                    auto result = iter->forward_label[1].get_f() +
                                  iter->y->backward_label[0].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    auto result = iter->forward_label[1].get_f() +
                                  iter->y->backward_label[1].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);

    return nb_removed_edges;
}

auto PricerSolverZddCycle::evaluate_nodes(double* pi) -> bool {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    auto reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (auto& it : table |
                        ranges::views::take(decision_diagram->topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::join) {
        auto* job = it.get_job();
        for (auto& iter : it.list) {
            auto w = iter->get_weight();

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (iter->forward_label[0].prev_job_forward() != job) {
                if (iter->y->backward_label[0].prev_job_backward() != job) {
                    auto result = iter->forward_label[0].get_f() +
                                  iter->y->backward_label[0].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    auto result = iter->forward_label[0].get_f() +
                                  iter->y->backward_label[1].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            } else {
                if (iter->y->backward_label[0].prev_job_backward() != job) {
                    auto result = iter->forward_label[1].get_f() +
                                  iter->y->backward_label[0].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    auto result = iter->forward_label[1].get_f() +
                                  iter->y->backward_label[1].get_f() -
                                  job->weighted_tardiness_start(w) +
                                  pi[job->job];
                    if (constLB - aux_nb_machines * reduced_cost - result >
                            UB + RC_FIXING &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);

    return nb_removed_edges;
}