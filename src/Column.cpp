#include "Column.h"                          // for ScheduleSet
#include <range/v3/algorithm/for_each.hpp>   // for for_each, for_each_fn
#include <range/v3/functional/identity.hpp>  // for identity
#include <utility>                           // for move
#include "Job.h"                             // for Job
#include "PricingSolution.hpp"               // for PricingSolution
#include "Solution.hpp"                      // for Machine

Column::Column(const Machine& m)
    : total_processing_time(m.total_processing_time),
      total_weighted_completion_time(m.total_weighted_tardiness),
      job_list(m.job_list) {}

Column::Column(PricingSolution&& pricing_solution)
    : total_processing_time(pricing_solution.C_max),
      total_weighted_completion_time(pricing_solution.cost),
      job_list(std::move(pricing_solution.jobs)) {}

// bool ScheduleSet::operator<(ScheduleSet const& other) {
// auto& tmp1 = job_list;
// auto& tmp2 = other.job_list;

// if (tmp1.size() != tmp2.size()) {
//     return tmp1.size() < tmp2.size();
// }

// for (auto i = 0UL; i < tmp1.size(); ++i) {
//     if (tmp1[i]->job != tmp2[i]->job) {
//         return tmp1[i]->job < tmp2[i]->job;
//     }
// }

// return job_list < other.job_list;
// }

// bool ScheduleSet::operator<(const ScheduleSet& other) {
//     return job_list < other.job_list;
// }

// bool ScheduleSet::operator==(ScheduleSet const& other) const {
// auto& tmp1 = job_list;
// auto& tmp2 = other.job_list;

// if (tmp1.size() != tmp2.size()) {
//     return false;
// }

// for (auto i = 0UL; i < tmp1.size(); ++i) {
//     if (tmp1[i]->job != tmp2[i]->job) {
//         return false;
//     }
// }

// return job_list == other.job_list;
// }

void Column::recalculate() {
    total_processing_time = 0;
    total_weighted_completion_time = 0;

    ranges::for_each(job_list, [&](Job* j) {
        total_processing_time += j->processing_time;
        total_weighted_completion_time +=
            j->weighted_tardiness(total_processing_time);
    });
}
