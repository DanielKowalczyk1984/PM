#include "PricerSolverZddForward.hpp"
#include "Instance.h"
#include "PricerSolverZdd.hpp"

/**
 *  zdd solver pricersolver for the flow formulation
 */
// PricerSolverSimple::PricerSolverSimple(GPtrArray*  _jobs,
//                                        int         _num_machines,
//                                        GPtrArray*  _ordered_jobs,
//                                        const char* _p_name,
//                                        double      _ub)
//     : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, _p_name, _ub) {
//     std::cout << "Constructing ZDD with Forward Simple evaluator" << '\n';
//     std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
//     std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
//     evaluator = ForwardZddSimpleDouble(convex_constr_id);
//     reversed_evaluator = BackwardZddSimpleDouble(convex_constr_id);
// }

PricerSolverSimple::PricerSolverSimple(const Instance& instance)
    : PricerSolverZdd(instance) {
    std::cout << "Constructing ZDD with Forward Simple evaluator" << '\n';
    std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
    evaluator = ForwardZddSimpleDouble(convex_constr_id);
    reversed_evaluator = BackwardZddSimpleDouble(convex_constr_id);
}

OptimalSolution<double> PricerSolverSimple::pricing_algorithm(double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

void PricerSolverSimple::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

void PricerSolverSimple::evaluate_nodes(double* pi, int UB, double LB) {
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
                        UB - 1 + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}

void PricerSolverSimple::evaluate_nodes(double* pi) {
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
                        UB - 1 + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}

// PricerSolverZddCycle::PricerSolverZddCycle(GPtrArray*  _jobs,
//                                            int         _num_machines,
//                                            GPtrArray*  _ordered_jobs,
//                                            const char* _p_name,
//                                            double      _ub)
//     : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, _p_name, _ub) {
//     std::cout << "Constructing ZDD with Forward ZddCycle evaluator" << '\n';
//     std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
//     std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
//     evaluator = ForwardZddCycleDouble(convex_constr_id);
//     reversed_evaluator = BackwardZddCycleDouble(convex_constr_id);
// }

PricerSolverZddCycle::PricerSolverZddCycle(const Instance& instance)
    : PricerSolverZdd(instance) {
    std::cout << "Constructing ZDD with Forward ZddCycle evaluator" << '\n';
    std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
    evaluator = ForwardZddCycleDouble(convex_constr_id);
    reversed_evaluator = BackwardZddCycleDouble(convex_constr_id);
}

OptimalSolution<double> PricerSolverZddCycle::pricing_algorithm(double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

void PricerSolverZddCycle::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

void PricerSolverZddCycle::evaluate_nodes(double* pi, int UB, double LB) {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto* job = it.get_job();
            for (auto& iter : it.list) {
                auto w = iter->get_weight();
                auto p = job->processing_time;

                auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
                if (iter->forward_label[0].prev_job_forward() != job) {
                    if (iter->y->backward_label[0].prev_job_backward() != job) {
                        auto result = iter->forward_label[0].get_f() +
                                      iter->y->backward_label[0].get_f() -
                                      value_Fj(w + p, job) + pi[job->job] +
                                      pi[convex_constr_id];
                        if (LB - aux_nb_machines * reduced_cost - result >
                                UB + RC_FIXING &&
                            (iter->calc_yes)) {
                            iter->calc_yes = false;
                            nb_removed_edges++;
                        }
                    } else {
                        auto result = iter->forward_label[0].get_f() +
                                      iter->y->backward_label[1].get_f() -
                                      value_Fj(w + p, job) + pi[job->job] +
                                      pi[convex_constr_id];
                        if (LB - aux_nb_machines * reduced_cost - result >
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
                                      value_Fj(w + p, job) + pi[job->job] +
                                      pi[convex_constr_id];
                        if (LB - aux_nb_machines * reduced_cost - result >
                                UB + RC_FIXING &&
                            (iter->calc_yes)) {
                            iter->calc_yes = false;
                            nb_removed_edges++;
                        }
                    } else {
                        auto result = iter->forward_label[1].get_f() +
                                      iter->y->backward_label[1].get_f() -
                                      value_Fj(w + p, job) + pi[job->job] +
                                      pi[convex_constr_id];
                        if (LB - aux_nb_machines * reduced_cost - result >
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
    }
}

void PricerSolverZddCycle::evaluate_nodes(double* pi) {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    auto reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto* job = it.get_job();
            for (auto& iter : it.list) {
                auto w = iter->get_weight();
                auto p = job->processing_time;

                auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
                if (iter->forward_label[0].prev_job_forward() != job) {
                    if (iter->y->backward_label[0].prev_job_backward() != job) {
                        auto result = iter->forward_label[0].get_f() +
                                      iter->y->backward_label[0].get_f() -
                                      value_Fj(w + p, job) + pi[job->job];
                        if (constLB - aux_nb_machines * reduced_cost - result >
                                UB + RC_FIXING &&
                            (iter->calc_yes)) {
                            iter->calc_yes = false;
                            nb_removed_edges++;
                        }
                    } else {
                        auto result = iter->forward_label[0].get_f() +
                                      iter->y->backward_label[1].get_f() -
                                      value_Fj(w + p, job) + pi[job->job];
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
                                      value_Fj(w + p, job) + pi[job->job];
                        if (constLB - aux_nb_machines * reduced_cost - result >
                                UB + RC_FIXING &&
                            (iter->calc_yes)) {
                            iter->calc_yes = false;
                            nb_removed_edges++;
                        }
                    } else {
                        auto result = iter->forward_label[1].get_f() +
                                      iter->y->backward_label[1].get_f() -
                                      value_Fj(w + p, job) + pi[job->job];
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
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);
}