#include "PricerSolverBddBackward.hpp"
#include <fmt/core.h>
#include <range/v3/all.hpp>
#include "Instance.h"
#include "PricerSolverBdd.hpp"
#include "util.h"

/**
 * backward bdd pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */

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

void PricerSolverBddBackwardSimple::evaluate_nodes(double* pi) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    auto reduced_cost = table.node(1).forward_label[0].get_f();
    auto removed_edges = false;
    auto nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    // for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
    //     for (auto& it : table[i]) {
    for (auto& it :
         table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
             ranges::views ::drop(1) | ranges::views::reverse |
             ranges::views::join) {
        auto& child = table.node(it[1]);
        auto  result = it.forward_label[0].get_f() +
                      child.backward_label[0].get_f() + it.reduced_cost[1];

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
    // }

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

void PricerSolverBddBackwardCycle::evaluate_nodes(double* pi) {
    auto& table = *(get_decision_diagram().getDiagram());
    compute_labels(pi);
    auto reduced_cost =
        table.node(get_decision_diagram().root()).backward_label[0].get_f();
    auto removed_edges = false;
    auto nb_removed_edges_evaluate = 0;

    /** check for each node the Lagrangian dual */
    // for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
    //     for (auto& it : table[i]) {
    for (auto& it :
         table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
             ranges::views ::drop(1) | ranges::views::reverse |
             ranges::views::join) {
        auto* job = it.get_job();
        auto  result{0.0};
        auto& child = table.node(it[1]);

        if (it.forward_label[0].prev_job_forward() != job &&
            child.backward_label[0].prev_job_backward() != job) {
            result = it.forward_label[0].get_f() +
                     child.backward_label[0].get_f() + it.reduced_cost[1];

        } else if (it.forward_label[0].prev_job_forward() == job &&
                   child.backward_label[0].prev_job_backward() != job) {
            result = it.forward_label[1].get_f() +
                     child.backward_label[0].get_f() + it.reduced_cost[1];
        } else if (it.forward_label[0].prev_job_forward() != job &&
                   child.backward_label[0].prev_job_backward() == job) {
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
    // }

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
