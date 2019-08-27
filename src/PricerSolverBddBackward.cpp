#include "PricerSolverBddBackward.hpp"

/**
 * backward bdd pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
    GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs)
    : PricerSolverBdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing BDD with Backward Simple evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = BackwardBddSimpleDouble(njobs);
    reversed_evaluator = ForwardBddSimpleDouble(njobs);
}

OptimalSolution<double> PricerSolverBddBackwardSimple::pricing_algorithm(
    double* _pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverBddBackwardSimple::compute_labels(double* _pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double* pi, int UB,
                                                   double LB) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            int    w = it.get_weight();
            Job*   job = it.get_job();
            double result = it.forward_label[0].get_f() +
                            it.child[1]->backward_label[0].get_f() -
                            value_Fj(w + job->processing_time, job) +
                            pi[job->job] + pi[njobs];
            auto result_no = it.forward_label[0].get_f() +
                             it.child[0]->backward_label[0].get_f() + pi[njobs];

            if (LB - (double)(num_machines - 1) * reduced_cost - result >
                    UB - 1 + 0.0001 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                nb_removed_edges++;
            }

            if (LB - (double)(num_machines - 1) * reduced_cost - result_no >
                    UB - 1 + 0.0001 &&
                (it.calc_no)) {
                it.calc_no = false;
                nb_removed_edges++;
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

/**
 * Simple backward bdd pricersolver for the flow formulation
 */
PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
    GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs)
    : PricerSolverBdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing BDD with Backward Cycle evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = BackwardBddCycleDouble(njobs);
    reversed_evaluator = ForwardBddCycleDouble(njobs);
}

OptimalSolution<double> PricerSolverBddBackwardCycle::pricing_algorithm(
    double* _pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(double* _pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi, int UB,
                                                  double LB) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            int    w = it.get_weight();
            Job*   job = it.get_job();
            double result;

            if (it.forward_label[0].get_previous_job() != job &&
                it.child[1]->backward_label[0].get_prev_job() != job) {
                result = it.forward_label[0].get_f() +
                         it.child[1]->backward_label[0].get_f() -
                         value_Fj(w + job->processing_time, job) +
                         pi[job->job] + pi[njobs];

            } else if (it.forward_label[0].get_previous_job() == job &&
                       it.child[1]->backward_label[0].get_prev_job() != job) {
                result = it.forward_label[1].get_f() +
                         it.child[1]->backward_label[0].get_f() -
                         value_Fj(w + job->processing_time, job) +
                         pi[job->job] + pi[njobs];
            } else if (it.forward_label[0].get_previous_job() != job &&
                       it.child[1]->backward_label[0].get_prev_job() == job) {
                result = it.forward_label[0].get_f() +
                         it.child[1]->backward_label[1].get_f() -
                         value_Fj(w + job->processing_time, job) +
                         pi[job->job] + pi[njobs];
            } else {
                result = it.forward_label[1].get_f() +
                         it.child[1]->backward_label[1].get_f() -
                         value_Fj(w + job->processing_time, job) +
                         pi[job->job] + pi[njobs];
            }

            if (LB - (double)(num_machines - 1) * reduced_cost - result >
                    UB + 0.0001 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                nb_removed_edges++;
            }

            auto max = std::numeric_limits<double>::min();

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    auto result_no = it.forward_label[i].get_f() +
                                     it.child[0]->backward_label[j].get_f() +
                                     pi[njobs];
                    if (max < result_no) {
                        max = result_no;
                    }
                }
            }

            if (LB - (double)(num_machines - 1) * reduced_cost - max >
                    UB - 1 + 0.00001 &&
                (it.calc_no)) {
                it.calc_no = false;
                nb_removed_edges++;
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}
