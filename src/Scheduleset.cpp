
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

bool ScheduleSet::operator<(ScheduleSet const& other) {
    auto& tmp1 = job_list;
    auto& tmp2 = other.job_list;

    if (tmp1.size() != tmp2.size()) {
        return tmp1.size() < tmp2.size();
    }

    for (auto i = 0UL; i < tmp1.size(); ++i) {
        if (tmp1[i]->job != tmp2[i]->job) {
            return tmp1[i]->job < tmp2[i]->job;
        }
    }

    return true;
}

bool ScheduleSet::operator==(ScheduleSet const& other) const {
    auto& tmp1 = job_list;
    auto& tmp2 = other.job_list;

    if (tmp1.size() != tmp2.size()) {
        return false;
    }

    for (auto i = 0UL; i < tmp1.size(); ++i) {
        if (tmp1[i]->job != tmp2[i]->job) {
            return false;
        }
    }

    return true;
}

void ScheduleSet::recalculate() {
    total_processing_time = 0;
    total_weighted_completion_time = 0;

    std::ranges::for_each(job_list, [&](Job* j) {
        total_processing_time += j->processing_time;
        total_weighted_completion_time += value_Fj(total_processing_time, j);
    });
}
