#ifndef PRICER_SOLVER_BDD_FORWARD_HPP
#define PRICER_SOLVER_BDD_FORWARD_HPP
#include <memory>                 // for make_unique, unique_ptr
#include "Instance.h"             // for Instance
#include "NodeBdd.hpp"            // for NodeBdd
#include "PricerEvaluateBdd.hpp"  // for BackwardBddFarkasDouble, BackwardBd...
#include "PricerSolverBdd.hpp"    // for PricerSolverBdd
#include "PricingSolution.hpp"    // for PricingSolution
struct PricerSolverBase;

class PricerSolverBddSimple : public PricerSolverBdd {
   private:
    ForwardBddSimpleDouble  evaluator;
    BackwardBddSimpleDouble reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddSimple(const Instance& instance);

    PricerSolverBddSimple(const PricerSolverBddSimple&) = default;
    PricerSolverBddSimple(PricerSolverBddSimple&&) = default;
    PricerSolverBddSimple& operator=(PricerSolverBddSimple&&) = default;
    PricerSolverBddSimple& operator=(const PricerSolverBddSimple&) = delete;
    ~PricerSolverBddSimple() override = default;

    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverBddSimple>(*this);
    };

    PricingSolution<double> pricing_algorithm(double* _pi) override;
    PricingSolution<double> pricing_algorithm(
        std::span<const double>& _pi) override;
    PricingSolution<double> farkas_pricing(double* _pi) override;
    PricingSolution<double> farkas_pricing(
        std::span<const double>& _pi) override;

    double evaluate_rc_arc(NodeBdd<>& n) override;
    void   compute_labels(double* _pi) override;
    void   compute_labels(std::span<const double>& _pi) override;
};

class PricerSolverBddCycle : public PricerSolverBdd {
   private:
    ForwardBddCycleDouble   evaluator;
    BackwardBddCycleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddCycle(const Instance& instance);

    PricerSolverBddCycle(const PricerSolverBddCycle&) = default;
    PricerSolverBddCycle(PricerSolverBddCycle&&) = default;
    PricerSolverBddCycle& operator=(PricerSolverBddCycle&&) = default;
    PricerSolverBddCycle& operator=(const PricerSolverBddCycle&) = delete;
    ~PricerSolverBddCycle() override = default;

    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverBddCycle>(*this);
    }

    PricingSolution<double> pricing_algorithm(double* _pi) override;
    PricingSolution<double> pricing_algorithm(
        std::span<const double>& _pi) override;
    PricingSolution<double> farkas_pricing(double* _pi) override;
    PricingSolution<double> farkas_pricing(
        std::span<const double>& _pi) override;

    double evaluate_rc_arc(NodeBdd<>& n) override;
    void   compute_labels(double* _pi) override;
    void   compute_labels(std::span<const double>& _pi) override;
};

#endif  // PRICER_SOLVER_BDD_FORWARD_HPP
