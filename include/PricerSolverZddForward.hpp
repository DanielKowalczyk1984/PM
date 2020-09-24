#ifndef PRICER_SOLVER_ZDD_FORWARD_HPP
#define PRICER_SOLVER_ZDD_FORWARD_HPP
#include "PricerEvaluateZdd.hpp"
#include "PricerSolverZdd.hpp"

class PricerSolverSimple : public PricerSolverZdd {
   private:
    ForwardZddSimpleDouble  evaluator;
    BackwardZddSimpleDouble reversed_evaluator;

   public:
    PricerSolverSimple(GPtrArray*  _jobs,
                       int         _num_machines,
                       GPtrArray*  _ordered_jobs,
                       const char* p_name,
                       double      _UB);
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};

class PricerSolverZddCycle : public PricerSolverZdd {
   private:
    ForwardZddCycleDouble  evaluator;
    BackwardZddCycleDouble reversed_evaluator;

   public:
    PricerSolverZddCycle(GPtrArray*  _jobs,
                         int         _num_machines,
                         GPtrArray*  _ordered_jobs,
                         const char* p_name,
                         double      _UB);
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};

#endif  // PRICER_SOLVER_ZDD_FORWARD_HPP
