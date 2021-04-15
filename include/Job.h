#ifndef JOB_H
#define JOB_H

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>

struct Job {
    int job{};
    int weight{};
    int processing_time{};
    int release_time{};
    int due_time{};

    Job() = default;
    Job(int p, int w, int d);
    Job(const Job&) = default;
    Job(Job&&) = default;
    Job& operator=(const Job&) = default;
    Job& operator=(Job&&) = default;
    ~Job() = default;

    int weighted_tardiness(int C);
    int weighted_tardiness_start(int S);
};

namespace std {
template <>
struct less<Job*> {
    constexpr bool operator()(auto const& lhs, auto const& rhs) {
        return (*lhs)->job < (*rhs)->job;  // or use boost::hash_combine
    }
};

template <>
struct equal_to<Job*> {
    constexpr bool operator()(const auto& lhs, const auto& rhs) {
        return lhs->job == rhs->job;
    }
};

}  // namespace std

int value_diff_Fij(int C, Job* i, Job* j);
int bool_diff_Fij(int, Job*, Job*);

#endif  // JOB_H
