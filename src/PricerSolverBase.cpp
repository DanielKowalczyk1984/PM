#include "PricerSolverBase.hpp"

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
      problem_name(p_name),
      env(new GRBEnv()),
      model(new GRBModel(*env)) { }

PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   GPtrArray* _ordered_jobs, const char* p_name)
    : jobs(_jobs),
      nb_jobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(_ordered_jobs),
      nb_layers(ordered_jobs->len),
      problem_name(p_name),
      env(new GRBEnv()),
      model(new GRBModel(*env)) {}

PricerSolverBase::~PricerSolverBase() {}

void PricerSolverBase::print_num_paths() {}
