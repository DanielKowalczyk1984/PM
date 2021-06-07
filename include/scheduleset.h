#ifndef SCHEDULESET_H
#define SCHEDULESET_H
#include <fmt/ostream.h>
#include <functional>
#include <memory>
#include <vector>
#include "ModelInterface.hpp"
#include "OptimalSolution.hpp"
#include "Solution.hpp"

struct ScheduleSet {
    int age{};
    int del{};
    int total_processing_time{};
    int total_weighted_completion_time{};

    std::vector<Job*> job_list{};

    ScheduleSet() = default;
    explicit ScheduleSet(const Machine&);
    explicit ScheduleSet(OptimalSolution<double>&& pricing_solution);
    ~ScheduleSet() = default;
    ScheduleSet(ScheduleSet&&) = default;
    ScheduleSet& operator=(ScheduleSet&&) = default;
    ScheduleSet(const ScheduleSet&) = default;
    ScheduleSet& operator=(const ScheduleSet&) = default;

    bool operator<(ScheduleSet const& other);
    bool operator==(ScheduleSet const& other) const;
    bool operator<=>(ScheduleSet const& other);

    friend std::ostream& operator<<(std::ostream& os, const ScheduleSet& set) {
        for (auto& it : set.job_list) {
            os << it->job << " ";
        }
        os << "C = " << set.total_processing_time
           << " cost = " << set.total_weighted_completion_time << '\n';
        return os;
    }

    void recalculate();
};

// namespace std {
// template <>
// struct less<std::shared_ptr<ScheduleSet>> {
//     constexpr bool operator()(auto const& lhs, auto const& rhs) {
//         return (*lhs) < (*rhs);  // or use boost::hash_combine
//     }
// };
// template <>
// struct equal_to<std::shared_ptr<ScheduleSet>> {
//     constexpr bool operator()(auto const& lhs, auto const& rhs) {
//         return (*lhs) == (*rhs);  // or use boost::hash_combine
//     }
// };

// }  // namespace std

#endif  // SCHEDULESET_H
