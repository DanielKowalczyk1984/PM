#ifndef SEPARATION_SOLVER_HPP
#define SEPARATION_SOLVER_HPP

#include "ModelInterface.hpp"
#include "PricerSolverBase.hpp"
// #include "wctprivate.h"

class SeperationSolver {
   public:
    explicit SeperationSolver(PricerSolverBase* _solver);
    SeperationSolver(SeperationSolver&&) = default;
    SeperationSolver(const SeperationSolver&) = default;
    SeperationSolver& operator=(SeperationSolver&&) = default;
    SeperationSolver& operator=(const SeperationSolver&) = default;
    ~SeperationSolver();

   private:
    PricerSolverBase*            solver;
    std::vector<ConstraintBase*> cuts;

    int max_nb_cuts;
};

#endif