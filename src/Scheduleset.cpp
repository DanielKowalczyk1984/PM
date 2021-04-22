
////////////////////////////////////////////////////////////////
//                                                            //
//  scheduleset.c                                                //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <bits/ranges_algo.h>
#include <scheduleset.h>

ScheduleSet::ScheduleSet(const Machine& m)
    : job_list(m.job_list),
      total_processing_time(m.completion_time),
      total_weighted_completion_time(m.total_weighted_tardiness) {}

ScheduleSet::ScheduleSet(OptimalSolution<double>&& pricing_solution)
    : total_processing_time(pricing_solution.C_max),
      total_weighted_completion_time(pricing_solution.cost),
      job_list(std::move(pricing_solution.jobs)) {}

bool ScheduleSet::operator<(ScheduleSet const& other) {
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

    return job_list < other.job_list;
}

// bool ScheduleSet::operator<(const ScheduleSet& other) {
//     return job_list < other.job_list;
// }

bool ScheduleSet::operator==(ScheduleSet const& other) const {
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

    return job_list == other.job_list;
}

void ScheduleSet::recalculate() {
    total_processing_time = 0;
    total_weighted_completion_time = 0;

    std::ranges::for_each(job_list, [&](Job* j) {
        total_processing_time += j->processing_time;
        total_weighted_completion_time +=
            j->weighted_tardiness(total_processing_time);
    });
}