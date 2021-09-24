#ifndef PRICER_SOLVER_ZDD_FORWARD_HPP
#define PRICER_SOLVER_ZDD_FORWARD_HPP
#include "Instance.h"             // for Instance
#include "PricerEvaluateZdd.hpp"  // for BackwardZddCycleDouble, BackwardZdd...
#include "PricerSolverZdd.hpp"    // for PricerSolverZdd
#include "PricingSolution.hpp"    // for PricingSolution
class PricerSolverSimple : public PricerSolverZdd {
   private:
    ForwardZddSimpleDouble  evaluator;
    BackwardZddSimpleDouble reversed_evaluator;

   public:
    explicit PricerSolverSimple(const Instance& instance);
    PricerSolverSimple(PricerSolverSimple&&) = default;
    PricerSolverSimple(const PricerSolverSimple&) = default;
    PricerSolverSimple& operator=(const PricerSolverSimple&) = delete;
    PricerSolverSimple& operator=(PricerSolverSimple&&) = delete;

    PricingSolution<double> pricing_algorithm(double* _pi) override;
    PricingSolution<double> pricing_algorithm(
        std::span<const double>& _pi) override;
    void compute_labels(double* _pi);
    void compute_labels(std::span<const double>& _pi);
    bool evaluate_nodes(double* pi) final;
    bool evaluate_nodes(std::span<const double>& pi) final;
    ~PricerSolverSimple() override = default;
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
    PricingSolution<double> pricing_algorithm(double* _pi) override;
    PricingSolution<double> pricing_algorithm(
        std::span<const double>& _pi) override;
    void compute_labels(double* _pi);
    void compute_labels(std::span<const double>& _pi);
    bool evaluate_nodes(double* pi) final;
    bool evaluate_nodes(std::span<const double>& pi) final;
    PricerSolverZddCycle(const PricerSolverZddCycle&) = default;
    PricerSolverZddCycle(PricerSolverZddCycle&&) = default;
    PricerSolverZddCycle& operator=(PricerSolverZddCycle&&) = default;
    PricerSolverZddCycle& operator=(const PricerSolverZddCycle&) = delete;
    ~PricerSolverZddCycle() override = default;
};

#endif  // PRICER_SOLVER_ZDD_FORWARD_HPP
