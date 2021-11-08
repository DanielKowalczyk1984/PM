// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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

    PricingSolution pricing_algorithm(double* _pi) override;
    PricingSolution pricing_algorithm(std::span<const double>& _pi) override;
    void            compute_labels(double* _pi);
    void            compute_labels(std::span<const double>& _pi);
    bool            evaluate_nodes(double* pi) final;
    bool            evaluate_nodes(std::span<const double>& pi) final;
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
    PricingSolution pricing_algorithm(double* _pi) override;
    PricingSolution pricing_algorithm(std::span<const double>& _pi) override;
    void            compute_labels(double* _pi);
    void            compute_labels(std::span<const double>& _pi);
    bool            evaluate_nodes(double* pi) final;
    bool            evaluate_nodes(std::span<const double>& pi) final;
    PricerSolverZddCycle(const PricerSolverZddCycle&) = default;
    PricerSolverZddCycle(PricerSolverZddCycle&&) = default;
    PricerSolverZddCycle& operator=(PricerSolverZddCycle&&) = default;
    PricerSolverZddCycle& operator=(const PricerSolverZddCycle&) = delete;
    ~PricerSolverZddCycle() override = default;
};

#endif  // PRICER_SOLVER_ZDD_FORWARD_HPP
