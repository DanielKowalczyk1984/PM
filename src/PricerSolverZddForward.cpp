#include "PricerSolverZddForward.hpp"
#include <fmt/core.h>                            // for print
#include <array>                                 // for array
// #include <ext/alloc_traits.h>                    // for __alloc_traits<>::va...
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

PricingSolution PricerSolverSimple::pricing_algorithm(double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

PricingSolution PricerSolverSimple::pricing_algorithm(
    std::span<const double>& _pi) {
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

bool PricerSolverSimple::evaluate_nodes(std::span<const double>& pi) {
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

bool PricerSolverSimple::evaluate_nodes(double* pi) {
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

PricingSolution PricerSolverZddCycle::pricing_algorithm(double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

PricingSolution PricerSolverZddCycle::pricing_algorithm(
    std::span<const double>& _pi) {
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

bool PricerSolverZddCycle::evaluate_nodes(std::span<const double>& pi) {
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

bool PricerSolverZddCycle::evaluate_nodes(double* pi) {
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