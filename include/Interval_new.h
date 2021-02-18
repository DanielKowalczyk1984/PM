#ifndef __INTERVAL_H__
#define __INTERVAL_H__
#include <array>
#include <memory>
#include <vector>
#include "job.h"

struct Interval {
    int a{};
    int b{};
    int begin{};
    int key{};

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