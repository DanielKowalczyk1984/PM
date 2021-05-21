#include "PricerSolverZddBackward.hpp"
#include <fmt/core.h>
#include <range/v3/all.hpp>
#include "Instance.h"

/**
 *  bdd solver pricersolver for the flow formulation
 */

PricerSolverZddBackwardSimple::PricerSolverZddBackwardSimple(
    const Instance& instance)
    : PricerSolverZdd(instance) {
    std::cout << "Constructing ZDD with Backward Simple evaluator" << '\n';
    std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
    evaluator = BackwardZddSimpleDouble(convex_constr_id);
    reversed_evaluator = ForwardZddSimpleDouble(convex_constr_id);
}

OptimalSolution<double> PricerSolverZddBackwardSimple::pricing_algorithm(
    double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverZddBackwardSimple::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverZddBackwardSimple::evaluate_nodes(double* pi) {
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
}

// PricerSolverZddBackwardCycle::PricerSolverZddBackwardCycle(
//     GPtrArray*  _jobs,
//     int         _num_machines,
//     GPtrArray*  _ordered_jobs,
//     const char* _p_name,
//     double      _ub)
//     : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, _p_name, _ub) {
//     std::cout << "Constructing ZDD with Backward ZddCycle evaluator" << '\n';
//     std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
//     std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
//     evaluator = BackwardZddCycleDouble(convex_constr_id);
//     reversed_evaluator = ForwardZddCycleDouble(convex_constr_id);
// }

PricerSolverZddBackwardCycle::PricerSolverZddBackwardCycle(
    const Instance& instance)
    : PricerSolverZdd(instance) {
    std::cout << "Constructing ZDD with Backward ZddCycle evaluator" << '\n';
    std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
    evaluator = BackwardZddCycleDouble(convex_constr_id);
    reversed_evaluator = ForwardZddCycleDouble(convex_constr_id);
}

OptimalSolution<double> PricerSolverZddBackwardCycle::pricing_algorithm(
    double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverZddBackwardCycle::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverZddBackwardCycle::evaluate_nodes(double* pi) {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
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
            if (iter->forward_label[0].prev_job_forward() != job &&
                iter->y->backward_label[0].prev_job_backward() != job) {
                double result = iter->forward_label[0].get_f() +
                                iter->y->backward_label[0].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (iter->forward_label[0].prev_job_forward() == job &&
                       iter->y->backward_label[0].prev_job_backward() != job) {
                double result = iter->forward_label[1].get_f() +
                                iter->y->backward_label[0].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (iter->forward_label[0].prev_job_forward() != job &&
                       iter->y->backward_label[0].prev_job_backward() == job) {
                double result = iter->forward_label[0].get_f() +
                                iter->y->backward_label[1].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            } else {
                double result = iter->forward_label[1].get_f() +
                                iter->y->backward_label[1].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}
