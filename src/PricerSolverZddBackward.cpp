#include "PricerSolverZddBackward.hpp"
#include <fmt/core.h>

/**
 *  bdd solver pricersolver for the flow formulation
 */
PricerSolverZddBackwardSimple::PricerSolverZddBackwardSimple(
    GPtrArray*  _jobs,
    int         _num_machines,
    GPtrArray*  _ordered_jobs,
    const char* p_name,
    double      _UB)
    : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, p_name, _UB) {
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

void PricerSolverZddBackwardSimple::evaluate_nodes(double* pi,
                                                   int     UB,
                                                   double  LB) {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();

    nb_removed_edges = 0;

    // /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (auto& iter : it.list) {
                int    w = iter->get_weight();
                Job*   job = it.get_job();
                double result = iter->forward_label[0].get_f() +
                                iter->y->backward_label[0].get_f() -
                                value_Fj(w + job->processing_time, job) +
                                pi[job->job] + pi[convex_constr_id];
                auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
                if (LB - aux_nb_machines * reduced_cost - result >
                        UB - 1 + 0.0001 &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}

void PricerSolverZddBackwardSimple::evaluate_nodes(double* pi) {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();

    nb_removed_edges = 0;

    // /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (auto& iter : it.list) {
                int    w = iter->get_weight();
                Job*   job = it.get_job();
                double result = iter->forward_label[0].get_f() +
                                iter->y->backward_label[0].get_f() -
                                value_Fj(w + job->processing_time, job) +
                                pi[job->job];
                auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB - 1 + 0.0001 &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}

PricerSolverZddBackwardCycle::PricerSolverZddBackwardCycle(
    GPtrArray*  _jobs,
    int         _num_machines,
    GPtrArray*  _ordered_jobs,
    const char* p_name,
    double      _UB)
    : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, p_name, _UB) {
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
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto* job = it.get_job();
            auto  p = job->processing_time;
            for (auto& iter : it.list) {
                auto w = iter->get_weight();

                auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
                if (iter->forward_label[0].get_previous_job() != job &&
                    iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[0].get_f() +
                                    iter->y->backward_label[0].get_f() -
                                    value_Fj(w + p, job) + pi[job->job];
                    if (constLB + aux_nb_machines * reduced_cost + result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() == job &&
                           iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[1].get_f() +
                                    iter->y->backward_label[0].get_f() -
                                    value_Fj(w + p, job) + pi[job->job];
                    if (constLB + aux_nb_machines * reduced_cost + result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() != job &&
                           iter->y->backward_label[0].get_prev_job() == job) {
                    double result = iter->forward_label[0].get_f() +
                                    iter->y->backward_label[1].get_f() -
                                    value_Fj(w + p, job) + pi[job->job];
                    if (constLB + aux_nb_machines * reduced_cost + result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    double result = iter->forward_label[1].get_f() +
                                    iter->y->backward_label[1].get_f() -
                                    value_Fj(w + p, job) + pi[job->job];
                    if (constLB + aux_nb_machines * reduced_cost + result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}

void PricerSolverZddBackwardCycle::evaluate_nodes(double* pi,
                                                  int     UB,
                                                  double  LB) {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto* job = it.get_job();
            auto  p = job->processing_time;
            for (auto& iter : it.list) {
                int w = iter->get_weight();

                auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
                if (iter->forward_label[0].get_previous_job() != job &&
                    iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[0].get_f() +
                                    iter->y->backward_label[0].get_f() -
                                    value_Fj(w + p, job) + pi[job->job] +
                                    pi[convex_constr_id];
                    if (LB - aux_nb_machines * reduced_cost - result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() == job &&
                           iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[1].get_f() +
                                    iter->y->backward_label[0].get_f() -
                                    value_Fj(w + p, job) + pi[job->job] +
                                    pi[convex_constr_id];
                    if (LB - aux_nb_machines * reduced_cost - result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() != job &&
                           iter->y->backward_label[0].get_prev_job() == job) {
                    double result = iter->forward_label[0].get_f() +
                                    iter->y->backward_label[1].get_f() -
                                    value_Fj(w + p, job) + pi[job->job] +
                                    pi[convex_constr_id];
                    if (LB - aux_nb_machines * reduced_cost - result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    double result = iter->forward_label[1].get_f() +
                                    iter->y->backward_label[1].get_f() -
                                    value_Fj(w + p, job) + pi[job->job] +
                                    pi[convex_constr_id];
                    if (LB - aux_nb_machines * reduced_cost - result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}
