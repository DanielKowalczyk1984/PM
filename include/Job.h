#ifndef JOB_H
#define JOB_H

#include <algorithm>
#include <cstddef>
struct Job {
    int job{};
    int weight{};
    int processing_time{};
    int release_time{};
    int due_time{};
    int num_layers{};

    Job() = default;
    Job(int p, int w, int d);
    Job(const Job&) = default;
    Job(Job&&) = default;
    Job& operator=(const Job&) = default;
    Job& operator=(Job&&) = default;
    ~Job() = default;
};

inline int value_Fj(int C, Job* j) {
    return j->weight * std::max(0, C - j->due_time);
}
int value_diff_Fij(int C, Job* i, Job* j);
int bool_diff_Fij(int, Job*, Job*);
int arctime_diff_Fij(int weight, Job* i, Job* j);

#endif  // JOB_H
