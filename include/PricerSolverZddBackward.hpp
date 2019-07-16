#ifndef PRICER_SOLVER_ZDD_BACKWARD_HPP
#define PRICER_SOLVER_ZDD_BACKWARD_HPP

#include "PricerSolverZdd.hpp"
#include "PricerEvaluateZdd.hpp"

class PricerSolverZddBackwardSimple : public PricerSolverZdd {
private:
    BackwardZddSimpleDouble evaluator;
    ForwardZddSimpleDouble reversed_evaluator;
public:
    PricerSolverZddBackwardSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverZddBackwardCycle : public PricerSolverZdd {
private:
    BackwardZddCycleDouble evaluator;
    ForwardZddCycleDouble reversed_evaluator;
public:
    PricerSolverZddBackwardCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;        
};

#endif // PRICER_SOLVER_ZDD_BACKWARD_HPP
