#ifndef SCHEDULESET_H
#define SCHEDULESET_H

#include <ostream>              // for operator<<, ostream, basic_ostream::o...
#include <vector>               // for vector
#include "Job.h"                // for Job
#include "PricingSolution.hpp"  // for PricingSolution
struct Machine;

struct Column {
    int    age{};
    bool   del{};
    int    total_processing_time{};
    double total_weighted_completion_time{};

    std::vector<Job*> job_list{};

    Column() = default;
    explicit Column(const Machine&);
    explicit Column(PricingSolution<double>&& pricing_solution);
    ~Column() = default;
    Column(Column&&) = default;
    Column& operator=(Column&&) = default;
    Column(const Column&) = default;
    Column& operator=(const Column&) = default;

    bool operator<(Column const& other);
    bool operator==(Column const& other) const;
    bool operator<=>(Column const& other);

    friend std::ostream& operator<<(std::ostream& os, const Column& set) {
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
