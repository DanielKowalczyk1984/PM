#ifndef OPTIMAL_SOLUTION_HPP
#define OPTIMAL_SOLUTION_HPP
#include <iostream>
#include <job.h>

template <typename T>
class OptimalSolution {
public:
    T                obj;
    int              cost;
    int              C_max;
    GPtrArray       *jobs;
    GPtrArray       *e_list;

    /** Default constructor */
    OptimalSolution() : obj(0), cost(0), C_max(0), jobs(g_ptr_array_new()), e_list(g_ptr_array_new())
    {
    }

    explicit OptimalSolution(T _obj) : obj(_obj), cost(0), C_max(0), jobs(g_ptr_array_new()), e_list(g_ptr_array_new())
    {
    }

    /** Copy constructor */
    OptimalSolution(const OptimalSolution& other) :
        obj(other.obj), cost(other.cost), C_max(other.C_max),
        jobs(g_ptr_array_sized_new(other.jobs->len)),
        e_list(g_ptr_array_sized_new(other.e_list->len))
    {
        for (unsigned i = 0; i < other.jobs->len; ++i) {
            g_ptr_array_add(jobs, g_ptr_array_index(other.jobs, i));
        }

        for (unsigned i = 0; i < other.edges->len; ++i) {
            g_ptr_array_add(jobs, g_ptr_array_index(other.edges, i));
        }
    }

    /** Move constructor */
    OptimalSolution(OptimalSolution&& other)
    noexcept :  /* noexcept needed to enable optimizations in containers */
        obj(other.obj), cost(other.cost), C_max(other.C_max), jobs(other.jobs), e_list(other.e_list)
    {
        other.jobs = nullptr;
        other.e_list = nullptr;
    }

    /** Copy assignment operator */
    OptimalSolution& operator= (const OptimalSolution& other)
    {
        OptimalSolution tmp(other);         // re-use copy-constructor
        *this = move(tmp); // re-use move-assignment
        return *this;
    }

    /** Move assignment operator */
    OptimalSolution& operator= (OptimalSolution&& other) noexcept
    {
        obj = other.obj;
        cost = other.cost;
        C_max = other.C_max;
        g_ptr_array_free(jobs, TRUE);
        jobs = other.jobs;
        other.jobs = nullptr;
        g_ptr_array_free(e_list, TRUE);
        e_list = other.e_list;
        other.e_list = nullptr;
        return *this;
    }

    /** Destructor */
    ~OptimalSolution() noexcept
    {
        if (jobs) {
            g_ptr_array_free(jobs, TRUE);
        }

        if(e_list) {
            g_ptr_array_free(e_list, TRUE);
        }
    }

    friend std::ostream& operator<<(std::ostream&              os,
                                    OptimalSolution<T> const& o)
    {
        os << "obj = " << o.obj << "," << std::endl
           << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;
        g_ptr_array_foreach(o.jobs, g_print_machine, NULL);
        return os;
    };

    inline void push_job_back(Job *_job, double _pi) {
        g_ptr_array_add(jobs, _job);
        C_max += _job->processing_time;
        cost += value_Fj(C_max, _job);
        obj += -_pi + value_Fj(C_max, _job);
    }

    inline void push_job_back_farkas(Job *_job, double _pi) {
        g_ptr_array_add(jobs, _job);
        C_max += _job->processing_time;
        cost += value_Fj(C_max, _job);
        obj += -_pi;
    }

    inline void push_job_back(Job *_job, int C, double _pi) {
        g_ptr_array_add(jobs, _job);
        cost += value_Fj(C + _job->processing_time, _job);
        obj += -_pi + value_Fj(C + _job->processing_time, _job);
    }

};


#endif // OPTIMAL_SOLUTION_HPP



