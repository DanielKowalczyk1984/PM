#ifndef OPTIMAL_SOLUTION_HPP
#define OPTIMAL_SOLUTION_HPP
#include <iostream>
#include <range/v3/algorithm/reverse.hpp>
#include <vector>
#include "Job.h"

class PricingSolution {
   public:
    double            obj{};
    int               cost{};
    int               C_max{};
    std::vector<Job*> jobs{};

    /** Default constructor */
    PricingSolution() = default;

    explicit PricingSolution(double _obj) : obj(_obj) {}

    /** Copy constructor */
    PricingSolution(const PricingSolution& other) = delete;

    /** Move constructor */
    PricingSolution(PricingSolution&& other) noexcept = default;

    /** Copy assignment operator */
    PricingSolution& operator=(const PricingSolution& other) = delete;

    /** Move assignment operator */
    PricingSolution& operator=(PricingSolution&& other) noexcept = default;

    /** Destructor */
    ~PricingSolution() = default;

    friend std::ostream& operator<<(std::ostream&          os,
                                    PricingSolution const& o) {
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

    void reverse_jobs() { ranges::reverse(jobs); }
};

#endif  // OPTIMAL_SOLUTION_HPP
