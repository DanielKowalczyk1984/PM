#include "PricerSolverBase.hpp"
#include <memory>
#include <set>
#include <vector>
#include "PricerConstruct.hpp"
#include "boost/graph/graphviz.hpp"

/**
 * PricerSolverBase default COnstructor
 **/
PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   const char* p_name)
    : jobs(_jobs),
      nb_jobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(nullptr),
      nb_layers(0),
      problem_name(p_name)
// size_graph(0), nb_removed_edges(0), nb_removed_nodes(0)
{
    // decision_diagram = nullptr;
}

PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   GPtrArray* _ordered_jobs, const char* p_name)
    : jobs(_jobs),
      nb_jobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(_ordered_jobs),
      nb_layers(ordered_jobs->len),
      problem_name(p_name) {}

PricerSolverBase::~PricerSolverBase() {}

void PricerSolverBase::print_num_paths() {}

void PricerSolverBase::calculate_edges(ScheduleSet* set) {
    return;
}