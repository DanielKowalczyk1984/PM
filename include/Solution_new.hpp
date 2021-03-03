#ifndef __SOLUTION_NEW_H__
#define __SOLUTION_NEW_H__
#include <memory>
#include <vector>
#include "Interval_new.h"
#include "job.h"

using VecJobPtr = std::vector<std::shared_ptr<Job>>;
using VecIntervalPtr = std::vector<std::shared_ptr<Interval>>;

template <typename IT>
void swap_ranges(IT start_a, IT end_a, IT start_b, IT end_b) {
    auto it = std::rotate(start_a, start_b, end_b);
    auto new_start_a = (end_a - start_a) + it;
    std::rotate(it, new_start_a, end_b);
}

struct Machine {
    VecJobPtr job_list{};

    int completion_time{};
    int total_weighted_tardiness{};

    bool updated{false};

    Machine() = default;
    Machine(Machine&&) = default;
    Machine& operator=(Machine&&) = default;
    Machine& operator=(const Machine&) = default;
    Machine(const Machine&) = default;
    ~Machine() = default;

    void add_job(std::shared_ptr<Job> job);
    void reset_machine(std::vector<int>& c);
};

struct Sol {
    std::vector<Machine> machines{};

    int nb_jobs{};
    int nb_machines{};

    std::vector<int> c{};
    std::vector<int> u{};

    int tw{};
    int off{};

    Sol() = default;
    Sol(int _nb_jobs, int _nb_machines, int _off)
        : machines(_nb_machines),
          nb_jobs(_nb_jobs),
          nb_machines(_nb_machines),
          c(_nb_jobs, -1),
          u(_nb_jobs, -1),
          off(_off) {}

    void construct_edd(const VecJobPtr& v);
    void construct_spt(const VecJobPtr& v);
    void construct_random_fisher_yates(const VecJobPtr& v);
    void construct_random_shuffle(const VecJobPtr& v);
    void canonical_order(const VecIntervalPtr& intervals);

    void update_insertion_move(int i, int j, int k, int l);
    void update_swap_move(int i, int j, int k, int l1, int l2);
    void update_insertion_move_inter(int i, int j, int k, int l, int m);
    void update_swap_move_inter(int i, int j, int k1, int k2, int l1, int l2);
    void print_solution();

   private:
    void add_job_front_machine(const std::shared_ptr<Job>& job);
    static constexpr auto cmp_machines_completion =
        [](const auto& lhs, const auto& rhs) -> bool {
        return lhs.completion_time > rhs.completion_time;
    };
    void calculate_partition(const VecIntervalPtr& v);
    // static constexpr auto func =
};

#endif  // __SOLUTION_NEW_H__