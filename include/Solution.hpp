#ifndef __SOLUTION_NEW_H__
#define __SOLUTION_NEW_H__

#include <cstddef>  // for size_t
#include <limits>   // for numeric_limits
#include <memory>   // for shared_ptr
#include <random>   // for mt19937
#include <vector>   // for vector

/** Forward declaration */
struct Instance;  // lines 60-60
struct Interval;  // lines 19-19
struct Job;       // lines 17-17

using VecJobPtr = std::vector<std::shared_ptr<Job>>;
using VecJobRawPtr = std::vector<Job*>;
using VecIntervalPtr = std::vector<std::shared_ptr<Interval>>;

template <typename IT>
void swap_ranges(IT start_a, IT end_a, IT start_b, IT end_b) {
    auto it = std::rotate(start_a, start_b, end_b);
    auto new_start_a = (end_a - start_a) + it;
    std::rotate(it, new_start_a, end_b);
}

struct Machine {
    VecJobRawPtr job_list{};

    int total_processing_time{};
    int total_weighted_tardiness{};

    bool updated{true};

    Machine() = default;
    Machine(Machine&&) = default;
    Machine& operator=(Machine&&) = default;
    Machine& operator=(const Machine&) = default;
    Machine(const Machine&) = default;
    ~Machine() = default;

    void add_job(Job* job);
    void reset_machine(std::vector<int>& c);
};

struct Sol {
    std::vector<Machine> machines{};

    size_t nb_jobs{};
    size_t nb_machines{};

    std::vector<int>    c{};
    std::vector<size_t> u{};

    int tw{std::numeric_limits<int>::max()};
    int off{};

    Sol() = default;
    explicit Sol(const Instance& instance);
    Sol(const Sol&) = default;
    Sol& operator=(const Sol&) = default;
    Sol(Sol&&) = default;
    Sol& operator=(Sol&&) = default;
    ~Sol() = default;

    void construct_edd(VecJobPtr& v);
    void construct_spt(const VecJobPtr& v);
    void construct_random_fisher_yates(const VecJobPtr& v);
    void construct_random_shuffle(const VecJobPtr& v);
    void canonical_order(const VecIntervalPtr& intervals);
    void perturb_swap_inter(size_t l1, size_t l2, std::mt19937& mt);

    void update_insertion_move(size_t i, size_t j, size_t k, size_t l);
    void update_swap_move(size_t i, size_t j, size_t k, size_t l1, size_t l2);
    void update_insertion_move_inter(size_t i,
                                     size_t j,
                                     size_t k,
                                     size_t l,
                                     size_t m);
    void update_swap_move_inter(size_t i,
                                size_t j,
                                size_t k1,
                                size_t k2,
                                size_t l1,
                                size_t l2);
    void print_solution() const;

   private:
    void add_job_front_machine(Job* job);
    void calculate_partition(const VecIntervalPtr& v);
};

#endif  // __SOLUTION_NEW_H__