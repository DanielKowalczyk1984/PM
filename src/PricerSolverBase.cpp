#include "PricerSolverBase.hpp"
#include <memory>
#include <set>
#include <vector>
#include "PricerConstruct.hpp"
#include "boost/graph/graphviz.hpp"

/**
 * PricerSolverBase default COnstructor
 **/
PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines)
    : jobs(_jobs),
      njobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(nullptr),
      nlayers(0)
// size_graph(0), nb_removed_edges(0), nb_removed_nodes(0)
{
    // decision_diagram = nullptr;
}

PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   GPtrArray* _ordered_jobs)
    : jobs(_jobs),
      njobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(_ordered_jobs),
      nlayers(ordered_jobs->len) {}

// size_graph(0),
// nb_removed_edges(0), nb_removed_nodes(0),
// env(new GRBEnv()),
// model(new GRBModel(*env))
//     {
//     /**
//      * Construction of decision diagram
//      */
//     PricerConstruct ps(ordered_jobs);
//     decision_diagram = std::unique_ptr<DdStructure<> >(new
//     DdStructure<>(ps)); size_graph = decision_diagram->size();
// }

PricerSolverBase::~PricerSolverBase() {}

/**
 * Some getters
 */
// void PricerSolverBase::iterate_zdd() {
//     DdStructure<Node<double>>::const_iterator it = decision_diagram->begin();

//     for (; it != decision_diagram->end(); ++it) {
//         std::set<int>::const_iterator i = (*it).begin();

//         for (; i != (*it).end(); ++i) {
//             std::cout << nlayers - *i << " ";
//         }

//         std::cout << '\n';
//     }
// }

void PricerSolverBase::print_num_paths() {
    // cout << "Number of paths: " <<
    // decision_diagram->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

// void PricerSolverBase::create_dot_zdd(const char *name) {
//     std::ofstream file;
//     file.open(name);
//     decision_diagram->dumpDot(file);
//     file.close();
// }

// void PricerSolverBase::print_number_nodes_edges() {
//     printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
//            nb_removed_nodes);
// }

// int PricerSolverBase::get_num_remove_nodes() {
//     return nb_removed_nodes;
// }

// int PricerSolverBase::get_num_remove_edges() {
//     return nb_removed_edges;
// }

// size_t PricerSolverBase::get_datasize() {
//     return decision_diagram->size();
// }

// size_t PricerSolverBase::get_size_graph() {
//     return decision_diagram->size();
// }

// int PricerSolverBase::get_num_layers() {
//     return decision_diagram->topLevel();
// }

// double PricerSolverBase::get_cost_edge(int idx) {
//     return 0.0;
// }

/**
 * Reduced cost Fixing
 */
// void PricerSolverBase::evaluate_nodes(double *pi, int UB, double LB) {
//     return;
// }

void PricerSolverBase::calculate_edges(ScheduleSet* set) {
    return;
}