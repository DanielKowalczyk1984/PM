#ifndef __INSTANCE_H__
#define __INSTANCE_H__
#include <fmt/core.h>
#include <filesystem>
#include <limits>
#include <memory>
#include <vector>

#include "Interval_new.h"
#include "job.h"
#include "parms.h"
#include "util.h"

namespace fs = std::filesystem;
namespace std {
template <>
struct default_delete<Job> {
    void operator()(Job* ptr) { CC_IFFREE(ptr, Job); }
};
}  // namespace std

struct Instance {
    fs::path path_to_instance{};

    Parms* parms{};

    std::vector<std::shared_ptr<Job>>      jobs;
    std::vector<std::shared_ptr<Interval>> intervals;

    std::vector<std::pair<std::shared_ptr<Job>, std::shared_ptr<Interval>>>
        vector_pair;

    /** Summary of jobs */
    int nb_jobs{};
    int nb_machines{};
    int p_sum{};
    int pmax{0};
    int pmin{std::numeric_limits<int>::max()};
    int dmax{0};
    int dmin{std::numeric_limits<int>::max()};
    int H_min{};
    int H_max{};
    int off{};

    Instance(const fs::path& _path, Parms* _parms);

    Instance() = default;

   private:
    void calculate_H_max_H_min();
    void find_division();
};

#endif  // __INSTANCE_H__