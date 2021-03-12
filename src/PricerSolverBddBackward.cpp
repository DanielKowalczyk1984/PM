#include "PricerSolverBddBackward.hpp"
#include <fmt/core.h>
#include <iostream>
#include "Instance.h"
#include "PricerSolverBdd.hpp"
#include "util.h"

/**
 * backward bdd pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
// PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
//     GPtrArray*  _jobs,
//     int         _num_machines,
//     GPtrArray*  _ordered_jobs,
//     const char* _p_name,
//     int         _hmax,
//     int*        _take_jobs,
//     double      _ub)
//     : PricerSolverBdd(_jobs,
//                       _num_machines,
//                       _ordered_jobs,
//                       _p_name,
//                       _hmax,
//                       _take_jobs,
//                       _ub) {
//     if (dbg_lvl() > 0) {
//         fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:",
//         ALIGN,
//                    "Backward Simple Evaluator");
//         fmt::print(R"({0: <{1}}{2}
// )",
//                    "Number of vertices BDD", ALIGN, get_nb_vertices());
//         fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN,
//                    get_nb_edges());
//     }
//     evaluator.set_table(&(*(get_decision_diagram().getDiagram())));
// }

PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
    const Instance& instance)
    : PricerSolverBdd(instance) {
    if (dbg_lvl() > 0) {
        fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
                   "Backward Simple Evaluator");
        fmt::print(R"({0: <{1}}{2}
)",
                   "Number of vertices BDD", ALIGN, get_nb_vertices());
        fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN,
                   get_nb_edges());
    }
    evaluator.set_table(&(*(get_decision_diagram().getDiagram())));
}

OptimalSolution<double> PricerSolverBddBackwardSimple::pricing_algorithm(
    double* _pi) {
    evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_backward(evaluator);
}

OptimalSolution<double> PricerSolverBddBackwardSimple::farkas_pricing(
    double* _pi) {
    farkas_evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}
void PricerSolverBddBackwardSimple::compute_labels(double* _pi) {
    evaluator.set_pi(_pi);
    reversed_evaluator.set_pi(_pi);
    get_decision_diagram().compute_labels_backward(evaluator);
    get_decision_diagram().compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double* pi,
                                                   int     UB,
                                                   double  LB) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    std::span aux_pi{pi, reformulation_model.size()};
    double    reduced_cost =
        table.node(1).forward_label[0].get_f() + aux_pi[convex_constr_id];
    bool removed_edges = false;
    int  nb_edges_removed_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto&  child = table.node(it.branch[1]);
            double result = it.forward_label[0].get_f() +
                            child.backward_label[0].get_f() +
                            it.reduced_cost[1] + aux_pi[convex_constr_id];

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (LB + aux_nb_machines * reduced_cost + result > UB + RC_FIXING &&
                (it.calc[1])) {
                it.calc[1] = false;
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
                (it.calc[1])) {
                it.calc[1] = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        if (dbg_lvl() > 0) {
            fmt::print("Number of edges removed by evaluate nodes {0:<{1}}\n",
                       nb_removed_edges_evaluate, ALIGN_HALF);
            fmt::print("Total number of edges removed {0:<{1}}\n",
                       get_nb_removed_edges(), ALIGN_HALF);
            fmt::print("Number of edges {0:<{1}}\n", get_nb_edges(),
                       ALIGN_HALF);
        }
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
// PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
//     GPtrArray*  _jobs,
//     int         _num_machines,
//     GPtrArray*  _ordered_jobs,
//     const char* _p_name,
//     int         _hmax,
//     int*        _take_jobs,
//     double      _ub)
//     : PricerSolverBdd(_jobs,
//                       _num_machines,
//                       _ordered_jobs,
//                       _p_name,
//                       _hmax,
//                       _take_jobs,
//                       _ub) {
//     if (dbg_lvl() > 0) {
//         fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:",
//         ALIGN,
//                    "Backward Cycle Evaluator");
//         fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
//                    get_nb_vertices());
//         fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN,
//                    get_nb_edges());
//     }
//     evaluator.set_table(&(*(get_decision_diagram().getDiagram())));
// }

PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
    const Instance& instance)
    : PricerSolverBdd(instance) {
    if (dbg_lvl() > 0) {
        fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
                   "Backward Cycle Evaluator");
        fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
                   get_nb_vertices());
        fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN,
                   get_nb_edges());
    }
    evaluator.set_table(&(*(get_decision_diagram().getDiagram())));
}

OptimalSolution<double> PricerSolverBddBackwardCycle::pricing_algorithm(
    double* _pi) {
    evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_backward(evaluator);
}

OptimalSolution<double> PricerSolverBddBackwardCycle::farkas_pricing(
    double* _pi) {
    farkas_evaluator.set_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(double* _pi) {
    evaluator.set_pi(_pi);
    reversed_evaluator.set_pi(_pi);
    get_decision_diagram().compute_labels_backward(evaluator);
    get_decision_diagram().compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi,
                                                  int     UB,
                                                  double  LB) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(get_decision_diagram().root()).backward_label[0].get_f();
    bool removed_edges = false;
    int  nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            double result{};
            auto&  child = table.node(it.branch[1]);

            // if (it.forward_label[0].get_previous_job() != job &&
            //     it.child[1]->backward_label[0].get_prev_job() != job) {
            result = it.forward_label[0].get_f() +
                     child.backward_label[0].get_f() + it.reduced_cost[1];

            // } else if (it.forward_label[0].get_previous_job() == job &&
            //            it.child[1]->backward_label[0].get_prev_job() != job)
            //            {
            //     result = it.forward_label[1].get_f() +
            //              it.child[1]->backward_label[0].get_f() +
            //              it.reduced_cost[1];
            // } else if (it.forward_label[0].get_previous_job() != job &&
            //            it.child[1]->backward_label[0].get_prev_job() == job)
            //            {
            //     result = it.forward_label[0].get_f() +
            //              it.child[1]->backward_label[1].get_f() +
            //              it.reduced_cost[1];
            // } else {
            //     result = it.forward_label[1].get_f() +
            //              it.child[1]->backward_label[1].get_f() +
            //              it.reduced_cost[1];
            // }

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result >
                    UB + RC_FIXING &&
                (it.calc[1])) {
                it.calc[1] = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        if (dbg_lvl() > 0) {
            fmt::print("Number of edges removed by evaluate nodes {0: <{1}}\n",
                       nb_removed_edges_evaluate, ALIGN_HALF);
            fmt::print("Total number of edges removed {0: <{1}}\n",
                       get_nb_removed_edges(), ALIGN_HALF);
            fmt::print("Number of edges {0: <{1}}\n", get_nb_edges(),
                       ALIGN_HALF);
        }
        remove_layers();
        remove_edges();
        // init_table();
    }
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(get_decision_diagram().root()).backward_label[0].get_f();
    bool removed_edges = false;
    int  nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto*  job = it.get_job();
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
                    UB - 1.0 + RC_FIXING &&
                (it.calc[1])) {
                it.calc[1] = false;
                removed_edges = true;
                add_nb_removed_edges();
                nb_removed_edges_evaluate++;
            }
        }
    }

    if (removed_edges) {
        if (dbg_lvl() > 0) {
            fmt::print("{0: <{2}}{1}\n",
                       "Number of edges removed by evaluate "
                       "nodes",
                       nb_removed_edges_evaluate, ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                       get_nb_removed_edges(), ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Number of edges", get_nb_edges(),
                       ALIGN);
        }
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }
}
