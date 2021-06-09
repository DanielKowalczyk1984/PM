#include "PricerSolverBddForward.hpp"
#include <fmt/core.h>
#include <range/v3/all.hpp>
#include "Instance.h"
#include "PricerSolverBdd.hpp"

/**
 *  bdd solver pricersolver for the flow formulation
 */
PricerSolverBddSimple::PricerSolverBddSimple(const Instance& instance)
    : PricerSolverBdd(instance) {
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

double PricerSolverBddSimple::evaluate_rc_arc(NodeBdd<>& n) {
    auto& table = *(get_decision_diagram().getDiagram());

    auto& child = table.node(n[1]);
    return n.forward_label[0].get_f() + child.backward_label[0].get_f() +
           n.reduced_cost[1];
}

// bool PricerSolverBddSimple::evaluate_nodes(double* pi) {
//     auto& table = *(get_decision_diagram().getDiagram());
//     compute_labels(pi);
//     auto reduced_cost = table.node(1).forward_label[0].get_f();
//     auto removed_edges = false;
//     auto nb_removed_edges_evaluate = 0;

//     /** check for each node the Lagrangian dual */
//     for (auto& it :
//          table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
//              ranges::views ::drop(1) | ranges::views::reverse |
//              ranges::views::join) {
//         auto&  child = table.node(it[1]);
//         double result = it.forward_label[0].get_f() +
//                         child.backward_label[0].get_f() + it.reduced_cost[1];
//         auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
//         if (constLB + aux_nb_machines * reduced_cost + result >
//                 UB + RC_FIXING &&
//             (it.calc[1])) {
//             it.calc[1] = false;
//             add_nb_removed_edges();
//             removed_edges = true;
//             nb_removed_edges_evaluate++;
//         }
//     }

//     if (removed_edges) {
//         fmt::print("Number of edges removed by evaluate nodes {0:<{1}}\n",
//                    nb_removed_edges_evaluate, ALIGN_HALF);
//         fmt::print("Total number of edges removed {0:<{1}}\n",
//                    get_nb_removed_edges(), ALIGN_HALF);
//         fmt::print("Number of edges {0:<{1}}\n", get_nb_edges(), ALIGN_HALF);
//         remove_layers();
//         remove_edges();
//         bottum_up_filtering();
//         topdown_filtering();
//         cleanup_arcs();
//         construct_mipgraph();
//     }

//     return removed_edges;
// }

/**
 * bdd solver pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
PricerSolverBddCycle::PricerSolverBddCycle(const Instance& instance)
    : PricerSolverBdd(instance) {
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

double PricerSolverBddCycle::evaluate_rc_arc(NodeBdd<>& n) {
    auto& table = *(get_decision_diagram().getDiagram());
    auto* job = n.get_job();
    auto& child = table.node(n[1]);

    if (n.forward_label[0].prev_job_forward() != job &&
        child.backward_label[0].prev_job_backward() != job) {
        return n.forward_label[0].get_f() + child.backward_label[0].get_f() +
               n.reduced_cost[1];
    } else if (n.forward_label[0].prev_job_forward() == job &&
               child.backward_label[0].prev_job_backward() != job) {
        return n.forward_label[1].get_f() + child.backward_label[0].get_f() +
               n.reduced_cost[1];
    } else if (n.forward_label[0].prev_job_forward() != job &&
               child.backward_label[0].prev_job_backward() == job) {
        return n.forward_label[0].get_f() + child.backward_label[1].get_f() +
               n.reduced_cost[1];
    } else {
        return n.forward_label[1].get_f() + child.backward_label[1].get_f() +
               n.reduced_cost[1];
    }
}

// bool PricerSolverBddCycle::evaluate_nodes(double* pi) {
//     auto& table = *(get_decision_diagram().getDiagram());
//     compute_labels(pi);
//     auto reduced_cost = table.node(1).forward_label[0].get_f();
//     auto removed_edges = false;
//     auto nb_removed_edges_evaluate = 0;

//     /** check for each node the Lagrangian dual */
//     for (auto& it :
//          table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
//              ranges::views ::drop(1) | ranges::views::reverse |
//              ranges::views::join) {
//         auto* job = it.get_job();
//         auto  result{0.0};
//         auto& child = table.node(it[1]);

//         if (it.forward_label[0].prev_job_forward() != job &&
//             child.backward_label[0].prev_job_backward() != job) {
//             result = it.forward_label[0].get_f() +
//                      child.backward_label[0].get_f() + it.reduced_cost[1];
//         } else if (it.forward_label[0].prev_job_forward() == job &&
//                    child.backward_label[0].prev_job_backward() != job) {
//             result = it.forward_label[1].get_f() +
//                      child.backward_label[0].get_f() + it.reduced_cost[1];
//         } else if (it.forward_label[0].prev_job_forward() != job &&
//                    child.backward_label[0].prev_job_backward() == job) {
//             result = it.forward_label[0].get_f() +
//                      child.backward_label[1].get_f() + it.reduced_cost[1];
//         } else {
//             result = it.forward_label[1].get_f() +
//                      child.backward_label[1].get_f() + it.reduced_cost[1];
//         }

//         auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
//         if (constLB + aux_nb_machines * reduced_cost + result >
//                 UB + RC_FIXING &&
//             (it.calc[1])) {
//             it.calc[1] = false;
//             removed_edges = true;
//             add_nb_removed_edges();
//             nb_removed_edges_evaluate++;
//         }
//     }

//     if (removed_edges) {
//         fmt::print("Number of edges removed by evaluate nodes {0: <{1}}\n",
//                    nb_removed_edges_evaluate, ALIGN_HALF);
//         fmt::print("Total number of edges removed {0: <{1}}\n",
//                    get_nb_removed_edges(), ALIGN_HALF);
//         fmt::print("Number of edges {0: <{1}}\n", get_nb_edges(),
//         ALIGN_HALF); remove_layers(); remove_edges(); bottum_up_filtering();
//         topdown_filtering();
//         cleanup_arcs();
//         construct_mipgraph();
//     }

//     return removed_edges;
// }
