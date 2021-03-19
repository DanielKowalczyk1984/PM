#ifndef __INSTANCE_H__
#define __INSTANCE_H__
#include <bits/c++config.h>
#include <fmt/core.h>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <memory>
#include <vector>

#include "Interval.h"
#include "Job.h"
#include "Parms.h"

namespace fs = std::filesystem;

struct Instance {
    fs::path                               path_to_instance{};
    std::vector<std::shared_ptr<Job>>      jobs;
    std::vector<std::shared_ptr<Interval>> intervals;

    std::vector<std::pair<Job*, Interval*>> vector_pair;

    /** Summary of jobs */
    size_t nb_jobs{};
    size_t nb_machines{};
    int    p_sum{};
    int    pmax{};
    int    pmin{std::numeric_limits<int>::max()};
    int    dmax{};
    int    dmin{std::numeric_limits<int>::max()};
    int    H_min{};
    int    H_max{};
    int    off{};

    Instance(const Parms& _parms);

    Instance(const Instance&) = default;
    Instance& operator=(const Instance&) = delete;
    Instance& operator=(Instance&&) = delete;
    Instance(Instance&&) = default;
    ~Instance() = default;
    class InstanceException : public std::exception {
       public:
        InstanceException(const char* const msg = nullptr) : errmsg(msg) { }

        [[nodiscard]] const char* what() const noexcept override { return (errmsg); }

       private:
        const char* errmsg;
    };

   private:
    void calculate_H_max_H_min();
    void find_division();
};

#endif  // __INSTANCE_H__