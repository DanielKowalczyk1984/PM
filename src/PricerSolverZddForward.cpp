#include "PricerSolverZddForward.hpp"

/**
 *  zdd solver pricersolver for the flow formulation
 */
PricerSolverSimple::PricerSolverSimple(GPtrArray* _jobs, int _num_machines,
                                       GPtrArray* _ordered_jobs,const char* p_name)
    : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, p_name) {
    std::cout << "Constructing ZDD with Forward Simple evaluator" << '\n';
    std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
    evaluator = ForwardZddSimpleDouble(nb_jobs);
    reversed_evaluator = BackwardZddSimpleDouble(nb_jobs);
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
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
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
                                pi[job->job] + pi[nb_jobs];

                if (LB - (double)(num_machines - 1) * reduced_cost - result >
                        UB - 1 + 0.0001 &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

PricerSolverZddCycle::PricerSolverZddCycle(GPtrArray* _jobs, int _num_machines,
                                           GPtrArray* _ordered_jobs, const char* p_name)
    : PricerSolverZdd(_jobs, _num_machines, _ordered_jobs, p_name) {
    std::cout << "Constructing ZDD with Forward ZddCycle evaluator" << '\n';
    std::cout << "number vertices ZDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges ZDD = " << get_nb_edges() << '\n';
    evaluator = ForwardZddCycleDouble(nb_jobs);
    reversed_evaluator = BackwardZddCycleDouble(nb_jobs);
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
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (auto& iter : it.list) {
                int  w = iter->get_weight();
                Job* job = it.get_job();

                if (iter->forward_label[0].get_previous_job() != job &&
                    iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[0].get_f() +
                                    iter->y->backward_label[0].get_f() -
                                    value_Fj(w + job->processing_time, job) +
                                    pi[job->job] + pi[nb_jobs];
                    if (LB - (double)(num_machines - 1) * reduced_cost -
                                result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() == job &&
                           iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[1].get_f() +
                                    iter->y->backward_label[0].get_f() -
                                    value_Fj(w + job->processing_time, job) +
                                    pi[job->job] + pi[nb_jobs];
                    if (LB - (double)(num_machines - 1) * reduced_cost -
                                result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() != job &&
                           iter->y->backward_label[0].get_prev_job() == job) {
                    double result = iter->forward_label[0].get_f() +
                                    iter->y->backward_label[1].get_f() -
                                    value_Fj(w + job->processing_time, job) +
                                    pi[job->job] + pi[nb_jobs];
                    if (LB - (double)(num_machines - 1) * reduced_cost -
                                result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    double result = iter->forward_label[1].get_f() +
                                    iter->y->backward_label[1].get_f() -
                                    value_Fj(w + job->processing_time, job) +
                                    pi[job->job] + pi[nb_jobs];
                    if (LB - (double)(num_machines - 1) * reduced_cost -
                                result >
                            UB + 0.0001 &&
                        (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}