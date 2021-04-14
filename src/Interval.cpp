#include "Interval.h"
#include <bits/ranges_algo.h>
#include <cmath>
#include <memory>

Interval::Interval(int _a, int _b, const vector_ptr_jobs& _sigma)
    : a(_a),
      b(_b),
      sigma(_sigma) {
    auto cmp = compare_edd(a, b);
    std::ranges::sort(sigma, cmp);
}

IntervalPair::IntervalPair(int                         _a,
                           int                         _b,
                           const std::shared_ptr<Job>& j1,
                           const std::shared_ptr<Job>& j2,
                           Interval*                   interval)
    : left(_a),
      right(_b),
      jobs{j1, j2},
      I(interval) {}

int IntervalPair::operator()() {
    left = I->a;
    right = I->a + jobs[1]->processing_time;

    if (left > I->b - jobs[0]->processing_time) {
        return left;
    } else {
        if (value_diff_Fij(left, jobs[0].get(), jobs[1].get()) <= 0) {
            return left;
        }

        // for (int t = k - 1; t >= 0; t--) {
        //     tmp = (interval*)g_ptr_array_index(interval_array, t);
        //     pair->left = tmp->a + j->processing_time - i->processing_time;

        //     if (value_diff_Fij(pair->left, i, j) <= 0 && pair->left >= tmp->a
        //     &&
        //         pair->left <= tmp->b - i->processing_time) {
        //         break;
        //     }
        // }

        left = jobs[0]->due_time +
               static_cast<int>(
                   ceil(static_cast<double>(jobs[1]->weight *
                                            jobs[0]->processing_time) /
                        jobs[0]->weight)) -
               jobs[0]->processing_time;
        return left;
    }
}
