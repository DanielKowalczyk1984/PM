#ifndef __INSTANCE_H__
#define __INSTANCE_H__
#include <fmt/core.h>
#include <filesystem>
#include <limits>
#include <memory>
#include <vector>

#include "Interval.h"
#include "Job.h"
#include "Parms.h"

namespace fs = std::filesystem;

struct Instance {
    fs::path path_to_instance{};

    Parms* parms{};

    std::vector<std::shared_ptr<Job>>      jobs;
    std::vector<std::shared_ptr<Interval>> intervals;

    std::vector<std::pair<Job*, Interval*>> vector_pair;

    /** Summary of jobs */
    int nb_jobs{};
    int nb_machines{};
    int p_sum{};
    int pmax{};
    int pmin{std::numeric_limits<int>::max()};
    int dmax{};
    int dmin{std::numeric_limits<int>::max()};
    int H_min{};
    int H_max{};
    int off{};

    Instance(const fs::path& _path, Parms* _parms);
    Instance(Parms* _parms);

    Instance() = default;

    Instance(const Instance&) = default;
    Instance& operator=(const Instance&) = default;
    Instance& operator=(Instance&&) = default;
    Instance(Instance&&) = default;
    ~Instance() = default;

   private:
    void calculate_H_max_H_min();
    void find_division();
};

#endif  // __INSTANCE_H__