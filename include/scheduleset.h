#ifndef SCHEDULESET_H
#define SCHEDULESET_H
// #ifdef __cplusplus
// extern "C" {
// #endif
// #include <glib.h>
// #include <interval.h>
#include <memory>
#include <vector>
#include "Solution.hpp"

struct ScheduleSet {
    int age{};
    int del{};
    int total_processing_time{};
    int total_weighted_completion_time{};
    int id{-1};

    std::vector<Job*> job_list{};

    ScheduleSet() = default;
    explicit ScheduleSet(const Machine&);
    ~ScheduleSet() = default;
    ScheduleSet(ScheduleSet&&) = default;
    ScheduleSet& operator=(ScheduleSet&&) = default;
    ScheduleSet(const ScheduleSet&) = default;
    ScheduleSet& operator=(const ScheduleSet&) = default;

    bool operator<(ScheduleSet const& other);
    bool operator==(ScheduleSet const& other) const;

    void recalculate();
};

namespace std {
template <>
struct less<std::shared_ptr<ScheduleSet>> {
    constexpr bool operator()(auto const& lhs, auto const& rhs) {
        return (*lhs) < (*rhs);  // or use boost::hash_combine
    }
};

}  // namespace std

/*Sorting schedulesets*/
// int  scheduleset_less(ScheduleSet* c1, ScheduleSet* c2);
// gint g_scheduleset_less(gconstpointer a, gconstpointer b);
// void g_scheduleset_print(gpointer data, gpointer user_data);
// void g_sum_recalculate(gpointer data, gpointer user_data);

/** new approach for columns */
// void g_sum_processing_time(gpointer data, gpointer user_data);

#endif  // SCHEDULESET_H
