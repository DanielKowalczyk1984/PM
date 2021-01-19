#include "PricerSolverBddForward.hpp"
#include <span>
#include "fmt/core.h"

/**
 *  bdd solver pricersolver for the flow formulation
 */
PricerSolverBddSimple::PricerSolverBddSimple(GPtrArray*  _jobs,
                                             int         _num_machines,
                                             GPtrArray*  _ordered_jobs,
                                             const char* _p_name,
                                             int         _hmax,
                                             int*        _take_jobs,
                                             double      _ub)
    : PricerSolverBdd(_jobs,
                      _num_machines,
                      _ordered_jobs,
                      _p_name,
                      _hmax,
                      _take_jobs,
                      _ub) {
    fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
               "Forward Simple Evaluator");
    fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
               get_nb_vertices());
    fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN, get_nb_edges());
}

OptimalSolution<double> PricerSolverBddSimple::pricing_algorithm(double* _pi) {
    evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_forward(evaluator);
}

OptimalSolution<double> PricerSolverBddSimple::farkas_pricing(double* _pi) {
    farkas_evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddSimple::compute_labels(double* _pi) {
    evaluator.set_pi(_pi);
    reversed_evaluator.set_pi(_pi);
    get_decision_diagram().compute_labels_forward(evaluator);
    get_decision_diagram().compute_labels_backward(reversed_evaluator);
}

void PricerSolverBddSimple::evaluate_nodes(double* pi, int UB, double LB) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(1).forward_label[0].get_f() + pi[convex_constr_id];
    bool removed_edges = false;
    int  nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto&  child = table.node(it.branch[1]);
            double result = it.forward_label[0].get_f() +
                            child.backward_label[0].get_f() +
                            it.reduced_cost[1] + pi[convex_constr_id];
            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (LB + aux_nb_machines * reduced_cost + result > UB + RC_FIXING &&
                (it.calc_yes)) {
                it.calc_yes = false;
                add_nb_removed_edges();
                removed_edges = true;
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("Number of edges removed by evaluate nodes {{0}:<{1}}\n",
                   nb_removed_edges_evaluate, ALIGN_HALF);
        fmt::print("Total number of edges removed {{0}:<{1}}\n",
                   get_nb_removed_edges(), ALIGN_HALF);
        fmt::print("Number of edges {{0}:<{1}}\n", get_nb_edges(), ALIGN_HALF);
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}

void PricerSolverBddSimple::evaluate_nodes(double* pi) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    bool   removed_edges = false;
    int    nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto&  child = table.node(it.branch[1]);
            double result = it.forward_label[0].get_f() +
                            child.backward_label[0].get_f() +
                            it.reduced_cost[1];
            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result >
                    UB + RC_FIXING &&
                (it.calc_yes)) {
                it.calc_yes = false;
                add_nb_removed_edges();
                removed_edges = true;
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("Number of edges removed by evaluate nodes {0:<{1}}\n",
                   nb_removed_edges_evaluate, ALIGN_HALF);
        fmt::print("Total number of edges removed {0:<{1}}\n",
                   get_nb_removed_edges(), ALIGN_HALF);
        fmt::print("Number of edges {0:<{1}}\n", get_nb_edges(), ALIGN_HALF);
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}

/**
 * bdd solver pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
PricerSolverBddCycle::PricerSolverBddCycle(GPtrArray*  _jobs,
                                           int         _num_machines,
                                           GPtrArray*  _ordered_jobs,
                                           const char* _p_name,
                                           int         _hmax,
                                           int*        _take_jobs,
                                           double      _ub)
    : PricerSolverBdd(_jobs,
                      _num_machines,
                      _ordered_jobs,
                      _p_name,
                      _hmax,
                      _take_jobs,
                      _ub) {
    fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
               "Forward Cycle Evaluator");
    fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
               get_nb_vertices());
    fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN, get_nb_edges());
}

OptimalSolution<double> PricerSolverBddCycle::pricing_algorithm(double* _pi) {
    evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_forward(evaluator);
}

OptimalSolution<double> PricerSolverBddCycle::farkas_pricing(double* _pi) {
    farkas_evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddCycle::compute_labels(double* _pi) {
    evaluator.set_pi(_pi);
    reversed_evaluator.set_pi(_pi);
    get_decision_diagram().compute_labels_forward(evaluator);
    get_decision_diagram().compute_labels_backward(reversed_evaluator);
}

void PricerSolverBddCycle::evaluate_nodes(double* pi, int UB, double LB) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(1).forward_label[0].get_f() + pi[convex_constr_id];
    bool removed_edges = false;
    int  nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            Job*   job = it.get_job();
            double result{};
            auto&  child = table.node(it.branch[1]);

            if (it.forward_label[0].get_previous_job() != job &&
                child.backward_label[0].get_prev_job() != job) {
                result = it.forward_label[0].get_f() +
                         child.backward_label[0].get_f() + it.reduced_cost[1] +
                         pi[convex_constr_id];
            } else if (it.forward_label[0].get_previous_job() == job &&
                       child.backward_label[0].get_prev_job() != job) {
                result = it.forward_label[1].get_f() +
                         child.backward_label[0].get_f() + it.reduced_cost[1] +
                         pi[convex_constr_id];
            } else if (it.forward_label[0].get_previous_job() != job &&
                       child.backward_label[0].get_prev_job() == job) {
                result = it.forward_label[0].get_f() +
                         child.backward_label[1].get_f() + it.reduced_cost[1] +
                         pi[convex_constr_id];
            } else {
                result = it.forward_label[1].get_f() +
                         child.backward_label[1].get_f() + it.reduced_cost[1] +
                         pi[convex_constr_id];
            }

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (LB + aux_nb_machines * reduced_cost + result > UB + RC_FIXING &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }

            // auto max = std::numeric_limits<double>::min();

            // for (int i = 0; i < 2; i++) {
            //     for (int j = 0; j < 2; j++) {
            //         auto result_no = -it.forward_label[i].get_f() -
            //                          it.child[0]->backward_label[j].get_f() +
            //                          pi[nb_jobs];
            //         if (max < result_no) {
            //             max = result_no;
            //         }
            //     }
            // }

            // auto result_no = it.forward_label[0].get_f() +
            //                  it.child[0]->backward_label[0].get_f() +
            //                  pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }
            // result_no = it.forward_label[1].get_f() +
            //             it.child[0]->backward_label[0].get_f() + pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }
            // result_no = it.forward_label[0].get_f() +
            //             it.child[0]->backward_label[1].get_f() + pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }
            // result_no = it.forward_label[1].get_f() +
            //             it.child[0]->backward_label[1].get_f() + pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }

            // if (it.forward_label[0].get_previous_job() !=
            // it.child[0]->backward_label[0].get_prev_job()) {
            // auto result_no = it.forward_label[0].get_f() +
            //                  it.child[0]->backward_label[0].get_f() +
            //                  pi[nb_jobs];
            // auto min =(double)(num_machines - 1) * reduced_cost + result_no;
            // for (int i = 0; i < num_machines + 1; i++)
            // {
            //     if (min > )
            //     {
            //         /* code */
            //     }

            // }

            // if (LB - (double)(num_machines - 1) * reduced_cost - max >
            //         UB + 0.00001 &&
            //     (it.calc_no)) {
            //     it.calc_no = false;
            //     nb_removed_edges++;
            // }
        }
    }

    if (removed_edges) {
        fmt::print("Number of edges removed by evaluate nodes {0: <{1}}\n",
                   nb_removed_edges_evaluate, ALIGN_HALF);
        fmt::print("Total number of edges removed {0: <{1}}\n",
                   get_nb_removed_edges(), ALIGN_HALF);
        fmt::print("Number of edges {0: <{1}}\n", get_nb_edges(), ALIGN_HALF);
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}

void PricerSolverBddCycle::evaluate_nodes(double* pi) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label[0].get_f();
    bool   removed_edges = false;
    int    nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            Job*   job = it.get_job();
            double result{};
            auto&  child = table.node(it.branch[1]);

            if (it.forward_label[0].get_previous_job() != job &&
                child.backward_label[0].get_prev_job() != job) {
                result = it.forward_label[0].get_f() +
                         child.backward_label[0].get_f() + it.reduced_cost[1];
            } else if (it.forward_label[0].get_previous_job() == job &&
                       child.backward_label[0].get_prev_job() != job) {
                result = it.forward_label[1].get_f() +
                         child.backward_label[0].get_f() + it.reduced_cost[1];
            } else if (it.forward_label[0].get_previous_job() != job &&
                       child.backward_label[0].get_prev_job() == job) {
                result = it.forward_label[0].get_f() +
                         child.backward_label[1].get_f() + it.reduced_cost[1];
            } else {
                result = it.forward_label[1].get_f() +
                         child.backward_label[1].get_f() + it.reduced_cost[1];
            }

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result >
                    UB + RC_FIXING &&
                (it.calc_yes)) {
                it.calc_yes = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }

            // auto max = std::numeric_limits<double>::min();

            // for (int i = 0; i < 2; i++) {
            //     for (int j = 0; j < 2; j++) {
            //         auto result_no = -it.forward_label[i].get_f() -
            //                          it.child[0]->backward_label[j].get_f() +
            //                          pi[nb_jobs];
            //         if (max < result_no) {
            //             max = result_no;
            //         }
            //     }
            // }

            // auto result_no = it.forward_label[0].get_f() +
            //                  it.child[0]->backward_label[0].get_f() +
            //                  pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }
            // result_no = it.forward_label[1].get_f() +
            //             it.child[0]->backward_label[0].get_f() + pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }
            // result_no = it.forward_label[0].get_f() +
            //             it.child[0]->backward_label[1].get_f() + pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }
            // result_no = it.forward_label[1].get_f() +
            //             it.child[0]->backward_label[1].get_f() + pi[nb_jobs];
            // if (max < result_no) {
            //     max = result_no;
            // }

            // if (it.forward_label[0].get_previous_job() !=
            // it.child[0]->backward_label[0].get_prev_job()) {
            // auto result_no = it.forward_label[0].get_f() +
            //                  it.child[0]->backward_label[0].get_f() +
            //                  pi[nb_jobs];
            // auto min =(double)(it.child[1]->num_machines - 1) * reduced_cost
            // + result_no; for (int i = 0; i < num_machines + 1; i++)
            // {
            //     if (min > )
            //     {
            //         /* code */
            //     }

            // }

            // if (LB - (double)(num_machines - 1) * reduced_cost - max >
            //         UB + 0.00001 &&
            //     (it.calc_no)) {
            //     it.calc_no = false;
            //     nb_removed_edges++;
            // }
        }
    }

    if (removed_edges) {
        fmt::print("Number of edges removed by evaluate nodes {0: <{1}}\n",
                   nb_removed_edges_evaluate, ALIGN_HALF);
        fmt::print("Total number of edges removed {0: <{1}}\n",
                   get_nb_removed_edges(), ALIGN_HALF);
        fmt::print("Number of edges {0: <{1}}\n", get_nb_edges(), ALIGN_HALF);
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}