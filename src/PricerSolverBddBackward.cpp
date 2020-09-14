#include "PricerSolverBddBackward.hpp"
#include <fmt/core.h>
#include <iostream>

/**
 * backward bdd pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
    GPtrArray*  _jobs,
    int         _num_machines,
    GPtrArray*  _ordered_jobs,
    const char* p_name,
    int         _Hmax,
    int*        _take_jobs,
    double      _UB)
    : PricerSolverBdd(_jobs,
                      _num_machines,
                      _ordered_jobs,
                      p_name,
                      _Hmax,
                      _take_jobs,
                      _UB) {
    fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", 60,
               "Backward Simple Evaluator");
    fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", 60,
               get_nb_vertices());
    fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", 60, get_nb_edges());
}

OptimalSolution<double> PricerSolverBddBackwardSimple::pricing_algorithm(
    double* _pi) {
    evaluator.set_pi(_pi);
    return get_decision_diagram()->evaluate_backward(evaluator);
}

OptimalSolution<double> PricerSolverBddBackwardSimple::farkas_pricing(
    double* _pi) {
    farkas_evaluator.set_pi(_pi);
    return get_decision_diagram()->evaluate_backward(farkas_evaluator);
}
void PricerSolverBddBackwardSimple::compute_labels(double* _pi) {
    evaluator.set_pi(_pi);
    reversed_evaluator.set_pi(_pi);
    get_decision_diagram()->compute_labels_backward(evaluator);
    get_decision_diagram()->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double* pi,
                                                   int     UB,
                                                   double  LB) {
    NodeTableEntity<>& table =
        get_decision_diagram()->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost =
        table.node(1).forward_label[0].get_f() + pi[convex_constr_id];
    bool removed_edges = false;
    int  nb_edges_removed_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram()->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            double result = it.forward_label[0].get_f() +
                            it.child[1]->backward_label[0].get_f() +
                            it.reduced_cost[1] + pi[convex_constr_id];

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (LB + aux_nb_machines * reduced_cost + result > UB + 0.0001 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_edges_removed_evaluate++;
            }
        }
    }

    if (removed_edges) {
        // std::cout << "Number of edges removed by evaluate_nodes \t= "<<
        // nb_edges_removed_evaluate << "\n"; std::cout << "Total number of
        // edges removed \t\t\t= " << get_nb_removed_edges() << "\n";
        remove_layers();
        remove_edges();
        // init_table();
    }
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double* pi) {
    NodeTableEntity<>& table =
        get_decision_diagram()->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    bool   removed_edges = false;
    int    nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram()->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            double result = it.forward_label[0].get_f() +
                            it.child[1]->backward_label[0].get_f() +
                            it.reduced_cost[1];

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);

            if (constLB + aux_nb_machines * reduced_cost + result > UB + 1e-4 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("Number of edges removed by evaluate nodes {{0}:<{1}}\n",
                   nb_removed_edges_evaluate, 30);
        fmt::print("Total number of edges removed {{0}:<{1}}\n",
                   get_nb_removed_edges(), 30);
        fmt::print("Number of edges {{0}:<{1}}\n", get_nb_edges(), 30);
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}

/**
 * Simple backward bdd pricersolver for the flow formulation
 */
PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
    GPtrArray*  _jobs,
    int         _num_machines,
    GPtrArray*  _ordered_jobs,
    const char* p_name,
    int         _Hmax,
    int*        _take_jobs,
    double      _UB)
    : PricerSolverBdd(_jobs,
                      _num_machines,
                      _ordered_jobs,
                      p_name,
                      _Hmax,
                      _take_jobs,
                      _UB) {
    fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", 60,
               "Backward Cycle Evaluator");
    fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", 60,
               get_nb_vertices());
    fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", 60, get_nb_edges());
}

OptimalSolution<double> PricerSolverBddBackwardCycle::pricing_algorithm(
    double* _pi) {
    evaluator.set_pi(_pi);
    return get_decision_diagram()->evaluate_backward(evaluator);
}

OptimalSolution<double> PricerSolverBddBackwardCycle::farkas_pricing(
    double* _pi) {
    farkas_evaluator.set_pi(_pi);
    return get_decision_diagram()->evaluate_backward(farkas_evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(double* _pi) {
    evaluator.set_pi(_pi);
    reversed_evaluator.set_pi(_pi);
    get_decision_diagram()->compute_labels_backward(evaluator);
    get_decision_diagram()->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi,
                                                  int     UB,
                                                  double  LB) {
    NodeTableEntity<>& table =
        get_decision_diagram()->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost =
        table.node(get_decision_diagram()->root()).backward_label[0].get_f();
    bool removed_edges = false;
    int  nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram()->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            Job*   job = it.get_job();
            double result;

            // if (it.forward_label[0].get_previous_job() != job &&
            //     it.child[1]->backward_label[0].get_prev_job() != job) {
            result = it.forward_label[0].get_f() +
                     it.child[1]->backward_label[0].get_f() +
                     it.reduced_cost[1];

            // } else if (it.forward_label[0].get_previous_job() == job &&
            //            it.child[1]->backward_label[0].get_prev_job() != job)
            //            {
            //     result = it.forward_label[1].get_f() +
            //              it.child[1]->backward_label[0].get_f() +
            //              it.reduced_cost[1] + pi[nb_jobs];
            // } else if (it.forward_label[0].get_previous_job() != job &&
            //            it.child[1]->backward_label[0].get_prev_job() == job)
            //            {
            //     result = it.forward_label[0].get_f() +
            //              it.child[1]->backward_label[1].get_f() +
            //              it.reduced_cost[1] + pi[nb_jobs];
            // } else {
            //     result = it.forward_label[1].get_f() +
            //              it.child[1]->backward_label[1].get_f() +
            //              it.reduced_cost[1] + pi[nb_jobs];
            // }

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result > UB + 1e-4 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("Number of edges removed by evaluate nodes {{0}:<{1}}\n",
                   nb_removed_edges_evaluate, 30);
        fmt::print("Total number of edges removed {{0}:<{1}}\n",
                   get_nb_removed_edges(), 30);
        fmt::print("Number of edges {{0}:<{1}}\n", get_nb_edges(), 30);
        remove_layers();
        remove_edges();
        // init_table();
    }
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi) {
    NodeTableEntity<>& table =
        get_decision_diagram()->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost =
        table.node(get_decision_diagram()->root()).backward_label[0].get_f();
    bool removed_edges = false;
    int  nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram()->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            Job*   job = it.get_job();
            double result;

            // if (it.forward_label[0].get_previous_job() != job &&
            //     it.child[1]->backward_label[0].get_prev_job() != job) {
            result = it.forward_label[0].get_f() +
                     it.child[1]->backward_label[0].get_f() +
                     it.reduced_cost[1];

            // } else if (it.forward_label[0].get_previous_job() == job &&
            //            it.child[1]->backward_label[0].get_prev_job() != job)
            //            {
            //     result = it.forward_label[1].get_f() +
            //              it.child[1]->backward_label[0].get_f() +
            //              it.reduced_cost[1] + pi[nb_jobs];
            // } else if (it.forward_label[0].get_previous_job() != job &&
            //            it.child[1]->backward_label[0].get_prev_job() == job)
            //            {
            //     result = it.forward_label[0].get_f() +
            //              it.child[1]->backward_label[1].get_f() +
            //              it.reduced_cost[1] + pi[nb_jobs];
            // } else {
            //     result = it.forward_label[1].get_f() +
            //              it.child[1]->backward_label[1].get_f() +
            //              it.reduced_cost[1] + pi[nb_jobs];
            // }

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result > UB + 1e-2 &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("{0: <{2}}{1}\n",
                   "Number of edges removed by evaluate "
                   "nodes",
                   nb_removed_edges_evaluate, 60);
        fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                   get_nb_removed_edges(), 60);
        fmt::print("{0: <{2}}{1}\n", "Number of edges", get_nb_edges(), 60);
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}
