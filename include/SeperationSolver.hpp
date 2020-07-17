#ifndef SEPARATION_SOLVER_HPP
#define SEPARATION_SOLVER_HPP

#include "ModelInterface.hpp"
#include "PricerSolverBase.hpp"
#include "wctprivate.h"

class SeparationSolver {
    PricerSolverBase* pricer_solver;

   public:
    SeparationSolver(PricerSolverBase* _pricer_solver)
        : pricer_solver(_pricer_solver) {
        ReformulationModel* rmp = pricer_solver->get_reformulation_model();
    }
};

#endif