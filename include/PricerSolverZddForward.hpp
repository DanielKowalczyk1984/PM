#ifndef PRICER_SOLVER_ZDD_FORWARD_HPP
#define PRICER_SOLVER_ZDD_FORWARD_HPP
#include "Instance.h"             // for Instance
#include "OptimalSolution.hpp"    // for OptimalSolution
#include "PricerEvaluateZdd.hpp"  // for BackwardZddCycleDouble, BackwardZdd...
#include "PricerSolverZdd.hpp"    // for PricerSolverZdd
class PricerSolverSimple : public PricerSolverZdd {
   private:
    ForwardZddSimpleDouble  evaluator;
    BackwardZddSimpleDouble reversed_evaluator;

   public:
    // PricerSolverSimple(GPtrArray*  _jobs,
    //                    int         _num_machines,
    //                    GPtrArray*  _ordered_jobs,
    //                    const char* p_name,
    //                    double      _UB);
    explicit PricerSolverSimple(const Instance& instance);
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    void                    compute_labels(double* _pi);
    bool                    evaluate_nodes(double* pi) final;
    PricerSolverSimple(const PricerSolverSimple&) = default;
};

class PricerSolverZddCycle : public PricerSolverZdd {
   private:
    ForwardZddCycleDouble  evaluator;
    BackwardZddCycleDouble reversed_evaluator;

   public:
    // PricerSolverZddCycle(GPtrArray*  _jobs,
    //                      int         _num_machines,
    //                      GPtrArray*  _ordered_jobs,
    //                      const char* p_name,
    //                      double      _UB);
    explicit PricerSolverZddCycle(const Instance& instance);
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    void                    compute_labels(double* _pi);
    bool                    evaluate_nodes(double* pi) final;
    PricerSolverZddCycle(const PricerSolverZddCycle&) = default;
};

#endif  // PRICER_SOLVER_ZDD_FORWARD_HPP
