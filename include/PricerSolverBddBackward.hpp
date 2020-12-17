#include "PricerEvaluateBdd.hpp"
#include "PricerSolverBdd.hpp"

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
   private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    PricerSolverBddBackwardSimple(GPtrArray*  _jobs,
                                  int         _num_machines,
                                  GPtrArray*  _ordered_jobs,
                                  const char* p_name,
                                  int         _Hmax,
                                  int*        _take_jobs,
                                  double      _UB);
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
    PricerSolverBddBackwardSimple(const PricerSolverBddBackwardSimple& src,
                                  GPtrArray* _ordered_jobs)
        : PricerSolverBdd(src, _ordered_jobs),
          evaluator(src.evaluator),
          reversed_evaluator(src.reversed_evaluator),
          farkas_evaluator(src.farkas_evaluator){};
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
   private:
    BackwardBddCycleDouble  evaluator;
    ForwardBddCycleDouble   reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    PricerSolverBddBackwardCycle(GPtrArray*  _jobs,
                                 int         _num_machines,
                                 GPtrArray*  _ordered_jobs,
                                 const char* p_name,
                                 int         _Hmax,
                                 int*        _take_jobs,
                                 double      _UB);

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
    PricerSolverBddBackwardCycle(const PricerSolverBddBackwardCycle& src,
                                 GPtrArray* _ordered_jobs)
        : PricerSolverBdd(src, _ordered_jobs),
          evaluator(src.evaluator),
          reversed_evaluator(src.reversed_evaluator),
          farkas_evaluator(src.farkas_evaluator){};
};