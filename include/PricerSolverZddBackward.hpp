#ifndef PRICER_SOLVER_ZDD_BACKWARD_HPP
#define PRICER_SOLVER_ZDD_BACKWARD_HPP

#include "Instance.h"             // for Instance
#include "PricerEvaluateZdd.hpp"  // for BackwardZddCycleDouble, BackwardZdd...
#include "PricerSolverZdd.hpp"    // for PricerSolverZdd
#include "PricingSolution.hpp"    // for PricingSolution

class PricerSolverZddBackwardSimple : public PricerSolverZdd {
   private:
    BackwardZddSimpleDouble evaluator;
    ForwardZddSimpleDouble  reversed_evaluator;

   public:
    // PricerSolverZddBackwardSimple(GPtrArray*  _jobs,
    //                               int         _num_machines,
    //                               GPtrArray*  _ordered_jobs,
    //                               const char* p_name,
    //                               double      _UB);
    PricerSolverZddBackwardSimple(const Instance& instance);
    PricingSolution pricing_algorithm(double* _pi) override;
    PricingSolution pricing_algorithm(std::span<const double>& _pi) override;
    void            compute_labels(double* _pi);
    void            compute_labels(std::span<const double>& _pi);
    bool            evaluate_nodes(double* pi) final;
    bool            evaluate_nodes(std::span<const double>& pi) final;
    PricerSolverZddBackwardSimple(const PricerSolverZddBackwardSimple&) =
        default;
};

class PricerSolverZddBackwardCycle : public PricerSolverZdd {
   private:
    BackwardZddCycleDouble evaluator;
    ForwardZddCycleDouble  reversed_evaluator;

   public:
    // PricerSolverZddBackwardCycle(GPtrArray*  _jobs,
    //                              int         _num_machines,
    //                              GPtrArray*  _ordered_jobs,
    //                              const char* p_name,
    //                              double      _UB);
    PricerSolverZddBackwardCycle(const Instance& instance);

    PricerSolverZddBackwardCycle(const PricerSolverZddBackwardCycle&) = default;
    PricerSolverZddBackwardCycle(PricerSolverZddBackwardCycle&&) = default;
    PricerSolverZddBackwardCycle& operator=(PricerSolverZddBackwardCycle&&) =
        delete;
    PricerSolverZddBackwardCycle& operator=(
        const PricerSolverZddBackwardCycle&) = delete;
    PricingSolution pricing_algorithm(double* _pi) override;
    PricingSolution pricing_algorithm(std::span<const double>& _pi) override;
    void            compute_labels(double* _pi);
    void            compute_labels(std::span<const double>& _pi);
    bool            evaluate_nodes(double* pi) final;
    bool            evaluate_nodes(std::span<const double>& pi) final;
};

#endif  // PRICER_SOLVER_ZDD_BACKWARD_HPP
