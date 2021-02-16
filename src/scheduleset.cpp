
////////////////////////////////////////////////////////////////
//                                                            //
//  scheduleset.c                                                //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <defs.h>
#include <fmt/core.h>
#include <scheduleset.h>
#include <util.h>
#include <span>

void g_sum_processing_time(gpointer data, gpointer user_data) {
    Job*         j = static_cast<Job*>(data);
    ScheduleSet* set = static_cast<ScheduleSet*>(user_data);

    set->total_processing_time += j->processing_time;
    set->total_weighted_completion_time +=
        value_Fj(set->total_processing_time, j);
    g_ptr_array_add(set->job_list, j);
}

void g_sum_recalculate(gpointer data, gpointer user_data) {
    Job*         j = static_cast<Job*>(data);
    ScheduleSet* set = static_cast<ScheduleSet*>(user_data);

    set->total_processing_time += j->processing_time;
    set->total_weighted_completion_time +=
        value_Fj(set->total_processing_time, j);
}

void g_scheduleset_print(gpointer data, MAYBE_UNUSED gpointer user_data) {
    auto* tmp = static_cast<ScheduleSet*>(data);
    auto* tmp_a = tmp->job_list;
    fmt::print("Machine {}: ", tmp->id);

    g_ptr_array_foreach(tmp_a, g_print_machine, NULL);

    fmt::print("with C = {}, cost = {} and {} jobs\n",
               tmp->total_processing_time, tmp->total_weighted_completion_time,
               tmp_a->len);
}

ScheduleSet::ScheduleSet(GPtrArray* machine) {
    job_list = g_ptr_array_new();
    g_ptr_array_foreach(machine, g_sum_processing_time, this);
}

ScheduleSet::ScheduleSet(const ScheduleSet& other)
    : age(other.age),
      del(other.del),
      total_processing_time(other.total_processing_time),
      total_weighted_completion_time(other.total_weighted_completion_time),
      job_list(g_ptr_array_copy(other.job_list, nullptr, nullptr)),
      id(other.id) {}

ScheduleSet& ScheduleSet::operator=(const ScheduleSet& other) {
    if (this != &other) {
        if (job_list) {
            g_ptr_array_free(job_list, TRUE);
        }

        age = other.age;
        del = other.del;
        total_processing_time = other.total_processing_time;
        total_weighted_completion_time = other.total_weighted_completion_time;
        id = other.id;
        job_list = g_ptr_array_copy(other.job_list, nullptr, nullptr);
    }

    return *this;
}

bool ScheduleSet::operator<(ScheduleSet const& other) {
    auto* tmp1 = job_list;
    auto* tmp2 = other.job_list;

    if (tmp1->len != tmp2->len) {
        return tmp1->len < tmp2->len;
    }

    std::span span1{tmp1->pdata, tmp1->len};
    std::span span2{tmp2->pdata, tmp2->len};
    for (guint i = 0; i < tmp1->len; ++i) {
        auto* tmp_j1 = static_cast<Job*>(span1[i]);
        auto* tmp_j2 = static_cast<Job*>(span2[i]);
        if (tmp_j1->job != tmp_j2->job) {
            return tmp_j1->job < tmp_j2->job;
        }
    }

    return true;
}

bool ScheduleSet::operator==(ScheduleSet const& other) const {
    auto* tmp1 = job_list;
    auto* tmp2 = other.job_list;

    if (tmp1->len != tmp2->len) {
        return false;
    }

    std::span span1{tmp1->pdata, tmp1->len};
    std::span span2{tmp2->pdata, tmp2->len};
    for (guint i = 0; i < tmp1->len; ++i) {
        auto* tmp_j1 = static_cast<Job*>(span1[i]);
        auto* tmp_j2 = static_cast<Job*>(span2[i]);
        if (tmp_j1->job != tmp_j2->job) {
            return false;
        }
    }

    return true;
}

ScheduleSet::~ScheduleSet() {
    if (job_list) {
        g_ptr_array_free(job_list, TRUE);
    }
}

void ScheduleSet::recalculate() {
    total_processing_time = 0;
    total_weighted_completion_time = 0;

    g_ptr_array_foreach(job_list, g_sum_recalculate, this);
}
