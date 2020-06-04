#include "PricerEvaluateBdd.hpp"
#include "PricerSolverBdd.hpp"

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
   private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    PricerSolverBddBackwardSimple(GPtrArray* _jobs, int _num_machines,
                                  GPtrArray* _ordered_jobs, const char* p_name);
    PricerSolverBddBackwardSimple(GPtrArray* _jobs, int _num_machines,
                                  GPtrArray* _ordered_jobs, int* _take_jobs,
                                  int _Hmax, const char* _p_name);
    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
   private:
    BackwardBddCycleDouble evaluator;
    ForwardBddCycleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    PricerSolverBddBackwardCycle(GPtrArray* _jobs, int _num_machines,
                                 GPtrArray* _ordered_jobs, const char* p_name);
    PricerSolverBddBackwardCycle(GPtrArray* _jobs, int _num_machines,
                                  GPtrArray* _ordered_jobs, int* _take_jobs,
                                  int _Hmax, const char* _p_name);

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
};