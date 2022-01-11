// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
    PricerSolverZddBackwardSimple(const Instance& instance);
    PricerSolverZddBackwardSimple(const PricerSolverZddBackwardSimple&) =
        default;
    auto pricing_algorithm(double* _pi) -> PricingSolution override;
    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_nodes(double* pi) -> bool final;
    auto evaluate_nodes(std::span<const double>& pi) -> bool final;
    void compute_labels(double* _pi);
    void compute_labels(std::span<const double>& _pi);
};

class PricerSolverZddBackwardCycle : public PricerSolverZdd {
   private:
    BackwardZddCycleDouble evaluator;
    ForwardZddCycleDouble  reversed_evaluator;

   public:
    PricerSolverZddBackwardCycle(const Instance& instance);

    PricerSolverZddBackwardCycle(const PricerSolverZddBackwardCycle&) = default;
    PricerSolverZddBackwardCycle(PricerSolverZddBackwardCycle&&) = default;
    auto operator=(PricerSolverZddBackwardCycle&&)
        -> PricerSolverZddBackwardCycle& = delete;
    auto operator=(const PricerSolverZddBackwardCycle&)
        -> PricerSolverZddBackwardCycle& = delete;

    auto pricing_algorithm(double* _pi) -> PricingSolution override;
    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_nodes(double* pi) -> bool final;
    auto evaluate_nodes(std::span<const double>& pi) -> bool final;
    void compute_labels(double* _pi);
    void compute_labels(std::span<const double>& _pi);
};

#endif  // PRICER_SOLVER_ZDD_BACKWARD_HPP
