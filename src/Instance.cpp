#include "Instance.h"
#include <fmt/core.h>                            // for print
#include <algorithm>                             // for min, max, min_element
#include <cmath>                                 // for ceil, floor
#include <compare>                               // for operator<
#include <cstddef>                               // for size_t
#include <ext/alloc_traits.h>                    // for __alloc_traits<>::va...
#include <fstream>                               // for basic_istream::opera...
#include <functional>                            // for identity, __invoke
#include <memory>                                // for shared_ptr, __shared...
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=, operator+
#include <range/v3/numeric/accumulate.hpp>       // for accumulate, accumula...
#include <range/v3/view/reverse.hpp>             // for reverse_view, revers...
#include <range/v3/view/take.hpp>                // for take_view, take, tak...
#include <range/v3/view/view.hpp>                // for operator|, view_closure
#include <ranges>                                // for next
#include <string>                                // for getline, string
#include <tuple>                                 // for tie, operator<=>
#include <utility>                               // for move, pair, make_pair
#include <vector>                                // for vector, erase_if
#include "Interval.h"                            // for IntervalPair, Interval
#include "Job.h"                                 // for Job
#include "Parms.h"                               // for Parms

Instance::Instance(const Parms& _parms)
    : path_to_instance(_parms.jobfile),
      nb_machines(_parms.nb_machines) {
    std::ifstream in_file{path_to_instance};

    if (in_file.is_open()) {
        std::string str;
        if (getline(in_file, str)) {
            std::istringstream ss(str);
            ss >> nb_jobs;
        }

        int p{}, d{}, w{};
        int counter{};
        while (getline(in_file, str)) {
            std::istringstream ss(str);
            ss >> p >> d >> w;
            d = d / static_cast<int>(nb_machines);
            jobs.emplace_back(std::make_shared<Job>(p, w, d));
            auto* tmp = jobs.back().get();
            if (tmp->processing_time > tmp->due_time) {
                off += tmp->weight * (tmp->processing_time - tmp->due_time);
                tmp->due_time = tmp->processing_time;
            }

            p_sum += tmp->processing_time;
            pmin = std::min(pmin, tmp->processing_time);
            pmax = std::max(pmax, tmp->processing_time);
            dmin = std::min(dmin, tmp->due_time);
            dmax = std::max(dmax, tmp->due_time);

            tmp->job = counter++;
        }
        calculate_H_max_H_min();
        find_division();
    } else {
        throw InstanceException("Could not open file\n");
    }
}

void Instance::calculate_H_max_H_min() {
    auto temp = p_sum - pmax;
    auto temp_dbl = static_cast<double>(temp);
    temp_dbl = floor(temp_dbl / nb_machines);
    H_max = static_cast<int>(temp_dbl) + pmax;
    H_min = static_cast<int>(ceil(temp_dbl / nb_machines)) - pmax;

    std::ranges::sort(jobs, [](const auto& lhs, const auto& rhs) -> bool {
        return (lhs->processing_time < rhs->processing_time);
    });

    auto tmp = p_sum;
    H_min = p_sum;
    tmp = ranges::accumulate(
        jobs | ranges::views::reverse | ranges::views::take(nb_machines - 1),
        p_sum, std::minus<>{}, [](auto& it) { return it->processing_time; });

    H_min = static_cast<int>(ceil(tmp / nb_machines));
    fmt::print(
        R"(H_max = {}, H_min = {},  pmax = {}, pmin = {}, p_sum = {}, off = {}
)",
        H_max, H_min, pmax, pmin, p_sum, off);

    std::ranges::sort(jobs, [](const auto& x, const auto& y) -> bool {
        // if ((x->due_time > y->due_time)) {
        //     return (false);
        // } else if (x->due_time < y->due_time) {
        //     return (true);
        // } else if (x->processing_time > y->processing_time) {
        //     return (false);
        // } else if (x->processing_time < y->processing_time) {
        //     return (true);
        // } else if (x->weight < y->weight) {
        //     return (false);
        // } else if (x->weight > y->weight) {
        //     return (true);
        // } else if (x->job > y->job) {
        //     return (false);
        // } else {
        //     return (true);
        // }

        return std::tie(x->due_time, x->processing_time, x->weight, x->job) <
               std::tie(y->due_time, y->processing_time, y->weight, y->job);
    });

    auto index = 0UL;
    std::ranges::for_each(jobs, [&index](auto& it) { it->job = index++; });
}

void Instance::find_division() {
    int prev = 0;

    std::vector<std::shared_ptr<Interval>> tmp_ptr_vec{};
    for (auto i = 0UL; i < nb_jobs && prev < H_max; ++i) {
        auto tmp = std::min(H_max, jobs[i]->due_time);
        if (prev < tmp) {
            tmp_ptr_vec.push_back(std::make_shared<Interval>(prev, tmp, jobs));
            prev = jobs[i]->due_time;
        }
    }

    if (prev < H_max) {
        tmp_ptr_vec.push_back(std::make_shared<Interval>(prev, H_max, jobs));
    }

    for (auto& it : tmp_ptr_vec) {
        std::vector<IntervalPair> special_pairs{};
        for (auto i = 0U; i < it->sigma.size() - 1; i++) {
            auto job1 = it->sigma[i];
            for (size_t j = i + 1; j < it->sigma.size(); j++) {
                auto job2 = it->sigma[j];
                auto intervalpair =
                    IntervalPair(it->a, it->b, job1, job2, it.get());

                if (!(job1->due_time >= it->b) && !(job2->due_time >= it->b) &&
                    !((it->a + job2->processing_time >= it->b) ||
                      (intervalpair() <= it->a))) {
                    special_pairs.push_back(intervalpair);
                }
            }
        }

        if (!special_pairs.empty()) {
            std::vector<int> t{};
            t.push_back(it->a);

            while (!special_pairs.empty()) {
                auto min = std::min_element(
                    special_pairs.begin(), special_pairs.end(),
                    [](const auto& lhs, const auto& rhs) -> bool {
                        return (lhs.right < rhs.right);
                    });

                t.push_back(min->right);
                special_pairs.erase(min);

                std::erase_if(special_pairs, [&t](auto& it_tmp) -> bool {
                    auto tmp = t.back();
                    if ((tmp >= it_tmp.left && tmp <= it_tmp.right) ||
                        (tmp + it_tmp.jobs[1]->processing_time >=
                         it_tmp.I->b)) {
                        return true;
                    } else {
                        it_tmp.right = tmp + it_tmp.jobs[1]->processing_time;
                        return false;
                    }
                });
            }

            t.push_back(it->b);

            for (auto i = 1U; i < t.size(); i++) {
                intervals.push_back(
                    std::make_shared<Interval>(t[i - 1], t[i], jobs));
            }
        } else {
            intervals.push_back(std::make_shared<Interval>(it->a, it->b, jobs));
        }
    }

    for (auto& it : intervals) {
        for (auto& job : it->sigma) {
            if (job->processing_time <= it->b) {
                vector_pair.emplace_back(std::make_pair(job.get(), it.get()));
            }
        }
    }

    fmt::print(R"(The number of layers = {}
)",
               vector_pair.size());
}
