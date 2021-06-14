#include "SeperationSolver.hpp"
struct PricerSolverBase;

SeperationSolver::SeperationSolver(PricerSolverBase* _solver)
    : solver(_solver) {}
