#ifndef SEPARATION_SOLVER_HPP
#define SEPARATION_SOLVER_HPP

#include <vector>  // for vector
class ConstraintBase;
struct PricerSolverBase;
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