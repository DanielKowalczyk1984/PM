#ifndef PRICER_SOLVER_BDD_FORWARD_HPP
#define PRICER_SOLVER_BDD_FORWARD_HPP
#include "PricerEvaluateBdd.hpp"
#include "PricerSolverBdd.hpp"

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
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    compute_labels(double* _pi);
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
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};

#endif  // PRICER_SOLVER_BDD_FORWARD_HPP
