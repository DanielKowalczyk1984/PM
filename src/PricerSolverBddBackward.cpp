#include "PricerSolverBddBackward.hpp"
#include <iostream>

/**
 * backward bdd pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
    GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs,
    const char* p_name)
    : PricerSolverBdd(_jobs, _num_machines, _ordered_jobs, p_name) {
    std::cout << "Constructing BDD with Backward Simple evaluator:" << '\n';
    std::cout << "number vertices BDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges BDD = " << get_nb_edges() << '\n';
    evaluator = BackwardBddSimpleDouble(nb_jobs);
    reversed_evaluator = ForwardBddSimpleDouble(nb_jobs);
    farkas_evaluator = BackwardBddFarkasDouble(nb_jobs);
}

PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
    GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs, int* _take_jobs, int _Hmax,
    const char* p_name)
    : PricerSolverBdd(_jobs, _num_machines, _ordered_jobs,_take_jobs,_Hmax, p_name) {
    std::cout << "Constructing BDD with Backward Simple evaluator:" << '\n';
    std::cout << "number vertices BDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges BDD = " << get_nb_edges() << '\n';
    evaluator = BackwardBddSimpleDouble(nb_jobs);
    reversed_evaluator = ForwardBddSimpleDouble(nb_jobs);
}

OptimalSolution<double> PricerSolverBddBackwardSimple::pricing_algorithm(
    double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

OptimalSolution<double> PricerSolverBddBackwardSimple::farkas_pricing(
    double* _pi) {
    farkas_evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(farkas_evaluator);
}
void PricerSolverBddBackwardSimple::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double* pi, int UB,
                                                   double LB) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    bool removed_edges = false;
    int nb_edges_removed_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            int    w = it.get_weight();
            Job*   job = it.get_job();
            double result = it.forward_label[0].get_f() +
                            it.child[1]->backward_label[0].get_f() +
                            value_Fj(w + job->processing_time, job) -
                            pi[job->job] - pi[nb_jobs];

            if (LB + (double)(num_machines - 1) * reduced_cost + result >
                    UB + 0.0001 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                nb_removed_edges++;
                nb_edges_removed_evaluate++;
            }
        }
    }

    if(removed_edges) {
        std::cout << "Number of edges removed by evaluate_nodes = "<< nb_edges_removed_evaluate << "\n";
        std::cout << "Total number of edges removed " << nb_removed_edges << "\n";
        remove_layers();
        remove_edges();
        // init_table();
    }
}

/**
 * Simple backward bdd pricersolver for the flow formulation
 */
PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
    GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs,
    const char* p_name)
    : PricerSolverBdd(_jobs, _num_machines, _ordered_jobs, p_name) {
    std::cout << "Constructing BDD with Backward Cycle evaluator" << '\n';
    std::cout << "number vertices BDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges BDD = " << get_nb_edges() << '\n';
    evaluator = BackwardBddCycleDouble(nb_jobs);
    reversed_evaluator = ForwardBddCycleDouble(nb_jobs);
    farkas_evaluator = BackwardBddFarkasDouble(nb_jobs);
}

PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
    GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs, int* _take_jobs, int _Hmax,
    const char* p_name)
    : PricerSolverBdd(_jobs, _num_machines, _ordered_jobs,_take_jobs,_Hmax, p_name) {
    std::cout << "Constructing BDD with Backward Simple evaluator:" << '\n';
    std::cout << "number vertices BDD = " << get_nb_vertices() << '\n';
    std::cout << "number edges BDD = " << get_nb_edges() << '\n';
    evaluator = BackwardBddCycleDouble(nb_jobs);
    reversed_evaluator = ForwardBddCycleDouble(nb_jobs);
    farkas_evaluator = BackwardBddFarkasDouble(nb_jobs);
}

OptimalSolution<double> PricerSolverBddBackwardCycle::pricing_algorithm(
    double* _pi) {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

OptimalSolution<double> PricerSolverBddBackwardCycle::farkas_pricing(
    double* _pi) {
    farkas_evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(farkas_evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(double* _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi, int UB,
                                                  double LB) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    bool removed_edges = false;
    int nb_removed_edges_evaluate = 0;


    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            int    w = it.get_weight();
            Job*   job = it.get_job();
            double result;

            if (it.forward_label[0].get_previous_job() != job &&
                it.child[1]->backward_label[0].get_prev_job() != job) {
                result = it.forward_label[0].get_f() +
                         it.child[1]->backward_label[0].get_f() +
                         value_Fj(w + job->processing_time, job) -
                         pi[job->job] - pi[nb_jobs];

            } else if (it.forward_label[0].get_previous_job() == job &&
                       it.child[1]->backward_label[0].get_prev_job() != job) {
                result = it.forward_label[1].get_f() +
                         it.child[1]->backward_label[0].get_f() +
                         value_Fj(w + job->processing_time, job) -
                         pi[job->job] - pi[nb_jobs];
            } else if (it.forward_label[0].get_previous_job() != job &&
                       it.child[1]->backward_label[0].get_prev_job() == job) {
                result = it.forward_label[0].get_f() +
                         it.child[1]->backward_label[1].get_f() +
                         value_Fj(w + job->processing_time, job) -
                         pi[job->job] - pi[nb_jobs];
            } else {
                result = it.forward_label[1].get_f() +
                         it.child[1]->backward_label[1].get_f() +
                         value_Fj(w + job->processing_time, job) -
                         pi[job->job] - pi[nb_jobs];
            }

            if (LB + (double)(num_machines - 1) * reduced_cost + result >
                    UB + 0.0001 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                nb_removed_edges++;
                nb_removed_edges_evaluate++;
            }
        }
    }

    if(removed_edges) {
        std::cout << "Number of edges removed by evaluate_nodes = "<< nb_removed_edges_evaluate << "\n";
        std::cout << "Total number of edges removed " << nb_removed_edges << "\n";
        remove_layers();
        remove_edges();
        // init_table();
    }

}
