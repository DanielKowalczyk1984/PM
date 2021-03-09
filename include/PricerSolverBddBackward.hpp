#include <memory>
#include "Instance.h"
#include "PricerEvaluateBdd.hpp"
#include "PricerSolverBdd.hpp"
// #include "interval.h"

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
   private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    // PricerSolverBddBackwardSimple(GPtrArray*  _jobs,
    //                               int         _num_machines,
    //                               GPtrArray*  _ordered_jobs,
    //                               const char* p_name,
    //                               int         _Hmax,
    //                               int*        _take_jobs,
    //                               double      _UB);
    explicit PricerSolverBddBackwardSimple(const Instance& instance);

    PricerSolverBddBackwardSimple(const PricerSolverBddBackwardSimple&) =
        default;
    PricerSolverBddBackwardSimple(PricerSolverBddBackwardSimple&&) = default;
    PricerSolverBddBackwardSimple& operator=(PricerSolverBddBackwardSimple&&) =
        default;
    PricerSolverBddBackwardSimple& operator=(
        const PricerSolverBddBackwardSimple&) = delete;
    virtual ~PricerSolverBddBackwardSimple() override = default;

    std::unique_ptr<PricerSolverBase> clone() const override {
        // auto* tmp =
        //     g_ptr_array_copy(get_ordered_jobs(), g_copy_interval_pair, NULL);
        return std::make_unique<PricerSolverBddBackwardSimple>(*this);
    };

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;

    void compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
   private:
    BackwardBddCycleDouble  evaluator;
    ForwardBddCycleDouble   reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    // PricerSolverBddBackwardCycle(GPtrArray*  _jobs,
    //                              int         _num_machines,
    //                              GPtrArray*  _ordered_jobs,
    //                              const char* p_name,
    //                              int         _Hmax,
    //                              int*        _take_jobs,
    //                              double      _UB);

    explicit PricerSolverBddBackwardCycle(const Instance& instance);

    PricerSolverBddBackwardCycle& operator=(
        const PricerSolverBddBackwardCycle&) = delete;
    PricerSolverBddBackwardCycle& operator=(PricerSolverBddBackwardCycle&&) =
        default;
    PricerSolverBddBackwardCycle(PricerSolverBddBackwardCycle&&) = default;
    PricerSolverBddBackwardCycle(const PricerSolverBddBackwardCycle&) = default;
    virtual ~PricerSolverBddBackwardCycle() override = default;

    std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverBddBackwardCycle>(*this);
    };

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;

    void compute_labels(double* _pi);
    void evaluate_nodes(double* pi, int UB, double LB) override;
    void evaluate_nodes(double* pi) final;
};