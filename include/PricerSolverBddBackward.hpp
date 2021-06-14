// #include <memory>
// #include <vector>
// #include "Instance.h"
// #include "NodeBdd.hpp"
// #include "PricerEvaluateBdd.hpp"
// #include "PricerSolverBdd.hpp"
// #include "interval.h"
#include <memory>                 // for make_unique, unique_ptr
#include "Instance.h"             // for Instance
#include "NodeBdd.hpp"            // for NodeBdd
#include "OptimalSolution.hpp"    // for OptimalSolution
#include "PricerEvaluateBdd.hpp"  // for BackwardBddFarkasDouble, BackwardBd...
#include "PricerSolverBdd.hpp"    // for PricerSolverBdd
struct PricerSolverBase;

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
   private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddBackwardSimple(const Instance& instance);

    PricerSolverBddBackwardSimple(const PricerSolverBddBackwardSimple&) =
        default;
    PricerSolverBddBackwardSimple(PricerSolverBddBackwardSimple&&) = default;
    PricerSolverBddBackwardSimple& operator=(PricerSolverBddBackwardSimple&&) =
        default;
    PricerSolverBddBackwardSimple& operator=(
        const PricerSolverBddBackwardSimple&) = delete;
    ~PricerSolverBddBackwardSimple() override = default;

    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverBddBackwardSimple>(*this);
    };

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;

    double evaluate_rc_arc(NodeBdd<>& n) override;
    void   compute_labels(double* _pi) override;
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
   private:
    BackwardBddCycleDouble  evaluator;
    ForwardBddCycleDouble   reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddBackwardCycle(const Instance& instance);

    PricerSolverBddBackwardCycle& operator=(
        const PricerSolverBddBackwardCycle&) = delete;
    PricerSolverBddBackwardCycle& operator=(PricerSolverBddBackwardCycle&&) =
        default;
    PricerSolverBddBackwardCycle(PricerSolverBddBackwardCycle&&) = default;
    PricerSolverBddBackwardCycle(const PricerSolverBddBackwardCycle&) = default;
    ~PricerSolverBddBackwardCycle() override = default;

    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverBddBackwardCycle>(*this);
    };

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;

    double evaluate_rc_arc(NodeBdd<>& n) override;
    void   compute_labels(double* _pi) override;
};