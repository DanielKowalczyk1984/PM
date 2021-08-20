#include "PricerSolverBddBackward.hpp"
#include <fmt/core.h>                            // for print
#include <array>                                 // for array
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=
#include <range/v3/view/join.hpp>                // for join_view
#include <range/v3/view/view.hpp>                // for operator|
#include "Instance.h"                            // for Instance
#include "Label.hpp"                             // for Label
#include "NodeBdd.hpp"                           // for NodeBdd
#include "NodeBddStructure.hpp"                  // for DdStructure
#include "NodeBddTable.hpp"                      // for TableHandler, NodeTa...
#include "PricerSolverBase.hpp"                  // for PricerSolverBase::ALIGN
#include "PricerSolverBdd.hpp"                   // for PricerSolverBdd
#include "util.h"                                // for dbg_lvl

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

double PricerSolverBddBackwardSimple::evaluate_rc_arc(NodeBdd<>& node) {
    auto& table = *(get_decision_diagram().getDiagram());
    auto& child = table.node(node[1]);
    return node.forward_label[0].get_f() + child.backward_label[0].get_f() +
           node.get_reduced_cost()[1];
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

double PricerSolverBddBackwardCycle::evaluate_rc_arc(NodeBdd<>& s) {
    auto& table = *(get_decision_diagram().getDiagram());
    auto* job = s.get_job();
    auto& child = table.node(s[1]);

    if (s.forward_label[0].prev_job_forward() != job &&
        child.backward_label[0].prev_job_backward() != job) {
        return s.forward_label[0].get_f() + child.backward_label[0].get_f() +
               s.get_reduced_cost()[1];
    } else if (s.forward_label[0].prev_job_forward() == job &&
               child.backward_label[0].prev_job_backward() != job) {
        return s.forward_label[1].get_f() + child.backward_label[0].get_f() +
               s.get_reduced_cost()[1];
    } else if (s.forward_label[0].prev_job_forward() != job &&
               child.backward_label[0].prev_job_backward() == job) {
        return s.forward_label[0].get_f() + child.backward_label[1].get_f() +
               s.get_reduced_cost()[1];
    } else {
        return s.forward_label[1].get_f() + child.backward_label[1].get_f() +
               s.get_reduced_cost()[1];
    }
}
