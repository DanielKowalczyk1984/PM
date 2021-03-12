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
    OptimalSolution() = default;

    explicit OptimalSolution(T _obj) : obj(_obj) {}

    /** Copy constructor */
    OptimalSolution(const OptimalSolution& other) = delete;

    /** Move constructor */
    OptimalSolution(OptimalSolution&& other) = default;

    /** Copy assignment operator */
    OptimalSolution& operator=(const OptimalSolution& other) = delete;

    /** Move assignment operator */
    OptimalSolution& operator=(OptimalSolution&& other) = default;

    /** Destructor */
    ~OptimalSolution() = default;

    friend std::ostream& operator<<(std::ostream&             os,
                                    OptimalSolution<T> const& o) {
        os << "obj = " << o.obj << "," << std::endl
           << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;
        // g_ptr_array_foreach(o.jobs, g_print_machine, NULL);

        for (auto& it : o.jobs) {
            os << it->job << ' ';
        }
        //   [&os](const Job& tmp) { os << tmp.job << " "; });
        return os;
    };

    inline void push_job_back(Job* _job, double _pi) {
        jobs.push_back(_job);
        C_max += _job->processing_time;
        cost += value_Fj(C_max, _job);
        obj += _pi;
    }

    inline void push_job_back_farkas(Job* _job, double _pi) {
        jobs.push_back(_job);
        C_max += _job->processing_time;
        cost += value_Fj(C_max, _job);
        obj += -_pi;
    }

    inline void push_job_back(Job* _job, int C, double _pi) {
        jobs.push_back(_job);
        cost += value_Fj(C + _job->processing_time, _job);
        obj += _pi;
    }

    void reverse_jobs() { std::ranges::reverse(jobs); }
};

#endif  // OPTIMAL_SOLUTION_HPP
