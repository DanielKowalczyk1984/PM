#ifndef PRICER_SOLVER_BDD_FORWARD_HPP
#define PRICER_SOLVER_BDD_FORWARD_HPP
#include "PricerEvaluateBdd.hpp"
#include "PricerSolverBdd.hpp"
#include "interval.h"

class PricerSolverBddSimple : public PricerSolverBdd {
   private:
    ForwardBddSimpleDouble  evaluator;
    BackwardBddSimpleDouble reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    PricerSolverBddSimple(GPtrArray*  _jobs,
                          int         _num_machines,
                          GPtrArray*  _ordered_jobs,
                          const char* p_name,
                          int         _Hmax,
                          int*        _take_jobs,
                          double      _UB);
    PricerSolverBddSimple(const PricerSolverBddSimple& src,
                          GPtrArray*                   _ordered_jobs)
        : PricerSolverBdd(src, _ordered_jobs),
          evaluator(src.evaluator),
          reversed_evaluator(src.reversed_evaluator),
          farkas_evaluator(src.farkas_evaluator){};

    PricerSolverBddSimple(const PricerSolverBddSimple&) = default;
    PricerSolverBddSimple(PricerSolverBddSimple&&) = default;
    PricerSolverBddSimple& operator=(PricerSolverBddSimple&&) = default;
    PricerSolverBddSimple& operator=(const PricerSolverBddSimple&) = delete;
    virtual ~PricerSolverBddSimple() override = default;

    std::unique_ptr<PricerSolverBase> clone() override {
        auto* tmp =
            g_ptr_array_copy(get_ordered_jobs(), g_copy_interval_pair, NULL);
        return std::make_unique<PricerSolverBddSimple>(*this, tmp);
    };

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;

    void compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};

class PricerSolverBddCycle : public PricerSolverBdd {
   private:
    ForwardBddCycleDouble   evaluator;
    BackwardBddCycleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    PricerSolverBddCycle(GPtrArray*  _jobs,
                         int         _num_machines,
                         GPtrArray*  _ordered_jobs,
                         const char* p_name,
                         int         _Hmax,
                         int*        _take_jobs,
                         double      _UB);
    PricerSolverBddCycle(const PricerSolverBddCycle& src,
                         GPtrArray*                  _ordered_jobs)
        : PricerSolverBdd(src, _ordered_jobs),
          evaluator(src.evaluator),
          reversed_evaluator(src.reversed_evaluator),
          farkas_evaluator(src.farkas_evaluator){};

    PricerSolverBddCycle(const PricerSolverBddCycle&) = default;
    PricerSolverBddCycle(PricerSolverBddCycle&&) = default;
    PricerSolverBddCycle& operator=(PricerSolverBddCycle&&) = default;
    PricerSolverBddCycle& operator=(const PricerSolverBddCycle&) = delete;
    virtual ~PricerSolverBddCycle() override = default;

    std::unique_ptr<PricerSolverBase> clone() override {
        auto* tmp =
            g_ptr_array_copy(get_ordered_jobs(), g_copy_interval_pair, NULL);
        return std::make_unique<PricerSolverBddCycle>(*this, tmp);
    }

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;

    void compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};

#endif  // PRICER_SOLVER_BDD_FORWARD_HPP
