#ifndef PRICER_SOLVER_ZDD_FORWARD_HPP
#define PRICER_SOLVER_ZDD_FORWARD_HPP
#include "PricerSolverZdd.hpp"
#include "PricerEvaluateZdd.hpp"

class PricerSolverSimple : public PricerSolverZdd {
private:
    ForwardZddSimpleDouble evaluator;
    BackwardZddSimpleDouble reversed_evaluator;
public:
    PricerSolverSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverZddCycle : public PricerSolverZdd {
private:
    ForwardZddCycleDouble evaluator;
    BackwardZddCycleDouble reversed_evaluator;
public:
    PricerSolverZddCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;        
};


#endif // PRICER_SOLVER_ZDD_FORWARD_HPP
