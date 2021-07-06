#ifndef JOB_H
#define JOB_H

#include <cstddef>

struct Job {
    size_t job{};
    int    processing_time{};
    int    due_time{};
    int    weight{};
    int    release_time{};

    Job() = default;
    Job(int p, int w, int d);
    Job(const Job&) = default;
    Job(Job&&) = default;
    Job& operator=(const Job&) = default;
    Job& operator=(Job&&) = default;
    ~Job() = default;

    int weighted_tardiness(int C);
    int weighted_tardiness(size_t C);
    int weighted_tardiness_start(int S);
    int weighted_tardiness_start(size_t S);
    int weighted_tardiness_start(long S);
};

int value_diff_Fij(int C, Job* i, Job* j);
int bool_diff_Fij(int, Job*, Job*);

#endif  // JOB_H
