#include "Instance.h"
#include <fmt/core.h>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "Interval_new.h"
#include "interval.h"
#include "job.h"
#include "util.h"

Instance::Instance(const fs::path& _path, Parms* _parms)
    : path_to_instance(_path),
      parms(_parms),
      nb_machines(parms->nb_machines) {
    std::ifstream in_file{_path.c_str()};

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
            d = d / nb_machines;
            jobs.push_back(std::shared_ptr<Job>(job_alloc(&p, &w, &d),
                                                std::default_delete<Job>()));
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
    }
}

void Instance::calculate_H_max_H_min() {
    int    temp = 0;
    double temp_dbl = 0.0;

    temp = p_sum - pmax;
    temp_dbl = static_cast<double>(temp);
    temp_dbl = floor(temp_dbl / nb_machines);
    H_max = static_cast<int>(temp_dbl) + pmax;
    H_min = static_cast<int>(ceil(temp_dbl / nb_machines)) - pmax;

    std::ranges::sort(jobs, [](const auto& lhs, const auto& rhs) -> bool {
        return (lhs->processing_time < rhs->processing_time);
    });

    int    m = 0;
    int    i = nb_jobs - 1;
    double tmp = p_sum;
    H_min = p_sum;
    do {
        auto* job = jobs[i].get();
        tmp -= job->processing_time;
        m++;
        i--;

    } while (m < nb_machines - 1);

    H_min = static_cast<int>(ceil(tmp / nb_machines));
    fmt::print(
        R"(H_max = {}, H_min = {},  pmax = {}, pmin = {}, p_sum = {}, off = {}
)",
        H_max, H_min, pmax, pmin, p_sum, off);

    std::sort(jobs.begin(), jobs.end(),
              [](const auto& x, const auto& y) -> bool {
                  if (x->due_time > y->due_time) {
                      return (false);
                  } else if (x->due_time < y->due_time) {
                      return (true);
                  } else if (x->processing_time > y->processing_time) {
                      return (false);
                  } else if (x->processing_time < y->processing_time) {
                      return (true);
                  } else if (x->weight < y->weight) {
                      return (false);
                  } else if (x->weight > y->weight) {
                      return (true);
                  } else if (x->job > y->job) {
                      return (false);
                  } else {
                      return (true);
                  }
              });

    int index = 0;
    for (auto& it : jobs) {
        it->job = index++;
    }
}

void Instance::find_division() {
    int prev = 0;

    std::vector<std::shared_ptr<Interval>> tmp_ptr_vec{};
    for (int i = 0; i < nb_jobs && prev < H_max; ++i) {
        int tmp = std::min(H_max, jobs[i]->due_time);
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

                std::erase_if(special_pairs, [&t](auto& it) -> bool {
                    auto tmp = t.back();
                    if ((tmp >= it.left && tmp <= it.right) ||
                        (tmp + it.jobs[1]->processing_time >= it.I->b)) {
                        return true;
                    } else {
                        it.right = tmp + it.jobs[1]->processing_time;
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
        auto& tmp_sigma = it->sigma;
        for (auto& job : tmp_sigma) {
            if (job->processing_time <= it->b) {
                vector_pair.emplace_back(std::make_pair(job, it));
            }
        }
    }

    fmt::print(R"(The number of layers = {}
)",
               vector_pair.size());
}
