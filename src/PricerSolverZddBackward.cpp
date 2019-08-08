#include "PricerSolverZddBackward.hpp"

/**
 *  bdd solver pricersolver for the flow formulation
 */
PricerSolverZddBackwardSimple::PricerSolverZddBackwardSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    PricerSolverZdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing ZDD with Backward Simple evaluator" << '\n';
    std::cout << "size BDD = " << decision_diagram->size() << '\n';
    evaluator = BackwardZddSimpleDouble(njobs);
    reversed_evaluator = ForwardZddSimpleDouble(njobs);
}

OptimalSolution<double> PricerSolverZddBackwardSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverZddBackwardSimple::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverZddBackwardSimple::evaluate_nodes(double *pi, int UB, double LB) {
    NodeTableEntity<NodeZdd<>>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();

    nb_removed_edges = 0;

    // /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            for(auto &iter : it.list) {
                int w = iter->get_weight();
                Job *job = it.get_job();
                double result = iter->forward_label[0].get_f() + iter->y->backward_label[0].get_f() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];

                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}


PricerSolverZddBackwardCycle::PricerSolverZddBackwardCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) : 
    PricerSolverZdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing ZDD with Backward ZddCycle evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = BackwardZddCycleDouble(njobs);
    reversed_evaluator = ForwardZddCycleDouble(njobs);
}

OptimalSolution<double> PricerSolverZddBackwardCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverZddBackwardCycle::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverZddBackwardCycle::evaluate_nodes(double *pi, int UB, double LB) {
    NodeTableEntity<NodeZdd<>>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            for(auto &iter : it.list) {
                int w = iter->get_weight();
                Job *job = it.get_job();

                if(iter->forward_label[0].get_previous_job() != job && iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[0].get_f() + iter->y->backward_label[0].get_f() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                    if (LB - (double)(num_machines - 1)*reduced_cost - result > UB + 0.0001 && (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() == job && iter->y->backward_label[0].get_prev_job() != job) {
                    double result = iter->forward_label[1].get_f() + iter->y->backward_label[0].get_f() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                    if (LB - (double)(num_machines - 1)*reduced_cost - result > UB + 0.0001 && (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else if (iter->forward_label[0].get_previous_job() != job && iter->y->backward_label[0].get_prev_job() == job) {
                    double result = iter->forward_label[0].get_f() + iter->y->backward_label[1].get_f() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                    if (LB - (double)(num_machines - 1)*reduced_cost - result > UB + 0.0001 && (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                } else {
                    double result = iter->forward_label[1].get_f() + iter->y->backward_label[1].get_f() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                    if (LB - (double)(num_machines - 1)*reduced_cost - result > UB + 0.0001 && (iter->calc_yes)) {
                        iter->calc_yes = false;
                        nb_removed_edges++;
                    }
                }
            }
        }
    }
}