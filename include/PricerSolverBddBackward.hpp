#include "PricerSolverBdd.hpp"
#include "PricerEvaluateBdd.hpp"

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble reversed_evaluator;

public:
    PricerSolverBddBackwardSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
private:
    BackwardBddCycleDouble evaluator;
    ForwardBddCycleDouble reversed_evaluator;
public:
    PricerSolverBddBackwardCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);

    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};