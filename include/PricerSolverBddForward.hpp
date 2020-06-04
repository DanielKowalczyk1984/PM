#ifndef PRICER_SOLVER_BDD_FORWARD_HPP
#define PRICER_SOLVER_BDD_FORWARD_HPP
#include "PricerSolverBdd.hpp"
// #include "PricerEvaluate.hpp"
#include "PricerEvaluateBdd.hpp"

class PricerSolverBddSimple : public PricerSolverBdd {
private:
    ForwardBddSimpleDouble evaluator;
    BackwardBddSimpleDouble reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;
public:
    PricerSolverBddSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs, const char* p_name);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    OptimalSolution<double> farkas_pricing(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverBddCycle : public PricerSolverBdd {
private:
    ForwardBddCycleDouble evaluator;
    BackwardBddCycleDouble reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;
public:
    PricerSolverBddCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs, const char* p_name);
    OptimalSolution<double> pricing_algorithm(double *_pi) override;
    OptimalSolution<double> farkas_pricing(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;        
};

#endif // PRICER_SOLVER_BDD_FORWARD_HPP

