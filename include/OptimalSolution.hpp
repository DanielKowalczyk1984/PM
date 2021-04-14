#ifndef OPTIMAL_SOLUTION_HPP
#define OPTIMAL_SOLUTION_HPP
#include <bits/ranges_algo.h>
#include <fmt/core.h>
#include <iostream>
#include <span>
#include <vector>
#include "Job.h"

template <typename T = double>
class OptimalSolution {
   public:
    T                 obj{};
    int               cost{};
    int               C_max{};
    std::vector<Job*> jobs{};

    /** Default constructor */
    OptimalSolution<T>() = default;

    explicit OptimalSolution<T>(T _obj) : obj(_obj) {}

    /** Copy constructor */
    OptimalSolution<T>(const OptimalSolution<T>& other) = delete;

    /** Move constructor */
    OptimalSolution<T>(OptimalSolution<T>&& other) noexcept = default;

    /** Copy assignment operator */
    OptimalSolution<T>& operator=(const OptimalSolution<T>& other) = delete;

    /** Move assignment operator */
    OptimalSolution<T>& operator=(OptimalSolution<T>&& other) noexcept = default;

    /** Destructor */
    ~OptimalSolution() = default;

    friend std::ostream& operator<<(std::ostream&             os,
                                    OptimalSolution<T> const& o) {
        os << "obj = " << o.obj << "," << std::endl
           << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;

        for (auto& it : o.jobs) {
            os << it->job << ' ';
        }
        return os;
    };

    inline void push_job_back(Job* _job, double _pi) {
        jobs.emplace_back(_job);
        C_max += _job->processing_time;
        cost += _job->weighted_tardiness(C_max);
        obj += _pi;
    }

    inline void push_job_back(Job* _job, int C, double _pi) {
        jobs.emplace_back(_job);
        cost += _job->weighted_tardiness_start(C);
        obj += _pi;
    }

    void reverse_jobs() { std::ranges::reverse(jobs); }
};

#endif  // OPTIMAL_SOLUTION_HPP
