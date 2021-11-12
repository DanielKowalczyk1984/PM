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

#ifndef __INTERVAL_H__
#define __INTERVAL_H__
#include <array>    // for array
#include <compare>  // for operator>
#include <memory>   // for shared_ptr
#include <tuple>    // for operator<=>, tie
#include <vector>   // for vector
#include "Job.h"

struct compare_edd {
   private:
    int a;
    int b;

   public:
    compare_edd(int _a, int _b) : a(_a), b(_b){};

    bool operator()(const auto lhs, const auto rhs) {
        auto diff = b - a;
        auto w_lhs = (lhs->due_time >= b) ? 0.0
                                          : static_cast<double>(lhs->weight) /
                                                lhs->processing_time;
        auto w_rhs = (rhs->due_time >= b) ? 0.0
                                          : static_cast<double>(rhs->weight) /
                                                rhs->processing_time;

        if (lhs->processing_time >= diff) {
            if (rhs->processing_time < diff) {
                return true;
            } else {
                // if (w_lhs > w_rhs) {
                //     return true;
                // } else if (w_rhs > w_lhs) {
                //     return false;
                // } else if (lhs->processing_time > rhs->processing_time) {
                //     return true;
                // } else {
                //     return false;
                // }

                return std::tie(w_lhs, lhs->processing_time) >
                       std::tie(w_rhs, rhs->processing_time);
            }
        } else {
            if (rhs->processing_time >= diff) {
                return false;
            } else {
                // if (w_lhs > w_rhs) {
                //     return true;
                // } else if (w_rhs > w_lhs) {
                //     return false;
                // } else if (lhs->processing_time > rhs->processing_time) {
                //     return true;
                // } else {
                //     return false;
                // }
                return std::tie(w_lhs, lhs->processing_time) >
                       std::tie(w_rhs, rhs->processing_time);
            }
        }
    };
};

struct Interval {
    int a{};
    int b{};

    using vector_ptr_jobs = std::vector<std::shared_ptr<Job>>;

    vector_ptr_jobs sigma;

    Interval(int, int, const vector_ptr_jobs&);

    Interval() = default;
    ~Interval() = default;
    Interval(const Interval&) = default;
    Interval(Interval&&) = default;
    Interval& operator=(Interval&&) = default;
    Interval& operator=(const Interval&) = default;
};

struct IntervalPair {
    int left{};
    int right{};

    std::array<std::shared_ptr<Job>, 2> jobs{nullptr, nullptr};
    Interval*                           I{nullptr};
    IntervalPair(int                         _a,
                 int                         _b,
                 const std::shared_ptr<Job>& j1,
                 const std::shared_ptr<Job>& j2,
                 Interval*                   _I);

    int operator()();
};

#endif  // __INTERVAL_H__