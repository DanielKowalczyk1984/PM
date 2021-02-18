#include "Interval_new.h"
#include <fmt/core.h>
#include <memory>
#include "Instance.h"

struct compare_edd {
    int a;
    int b;
    compare_edd(int _a, int _b) : a(_a), b(_b){};

    bool operator()(const std::shared_ptr<Job> lhs,
                    const std::shared_ptr<Job> rhs) {
        int    diff = b - a;
        double w_lhs = (lhs->due_time >= b)
                           ? 0.0
                           : (double)lhs->weight / lhs->processing_time;
        double w_rhs = (rhs->due_time >= b)
                           ? 0.0
                           : (double)rhs->weight / rhs->processing_time;

        if (lhs->processing_time >= diff) {
            if (rhs->processing_time < diff) {
                return true;
            } else {
                if (w_lhs > w_rhs) {
                    return true;
                } else if (w_rhs > w_lhs) {
                    return false;
                } else if (lhs->processing_time > rhs->processing_time) {
                    return true;
                } else {
                    return false;
                }
            }
        } else {
            if (rhs->processing_time >= diff) {
                return false;
            } else {
                if (w_lhs > w_rhs) {
                    return true;
                } else if (w_rhs > w_lhs) {
                    return false;
                } else if (lhs->processing_time > rhs->processing_time) {
                    return true;
                } else {
                    return false;
                }
            }
        }
    }
};

Interval::Interval(int _a, int _b, const vector_ptr_jobs& _sigma)
    : a(_a),
      b(_b),
      begin(),
      key(),
      sigma(_sigma) {
    auto cmp = compare_edd(a, b);
    std::sort(sigma.begin(), sigma.end(), cmp);
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
