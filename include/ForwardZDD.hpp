#ifndef DURATION_ZDD_HPP
#define DURATION_ZDD_HPP
#include <algorithm>
#include <limits>
#include <vector>
#include "ForwardBDD.hpp"
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"
#include "ZddNode.hpp"

template <typename T = double>
class ForwardZddBase : public Eval<NodeZdd<T>, OptimalSolution<T>> {
   protected:
    T*  pi;
    int num_jobs;

   public:
    ForwardZddBase(T* _pi, int _num_jobs) : pi(_pi), num_jobs(_num_jobs) {}

    explicit ForwardZddBase(int _num_jobs)
        : Eval<NodeZdd<T>, OptimalSolution<T>>(),
          pi(nullptr),
          num_jobs(_num_jobs) {}

    ForwardZddBase() : Eval<NodeZdd<T>, OptimalSolution<T>>() {
        pi = nullptr;
        num_jobs = 0;
    }

    ~ForwardZddBase<T>() = default;

    // ForwardZddBase(const ForwardZddBase<T>& src) {
    //     pi = src.pi;
    //     num_jobs = src.num_jobs;
    // }

    void initialize_pi(T* _pi) { pi = _pi; }

    virtual void initializenode(NodeZdd<T>& n) const = 0;

    virtual void initializerootnode(NodeZdd<T>& n) const = 0;

    virtual void evalNode(NodeZdd<T>& n) const = 0;

    OptimalSolution<T> get_objective(NodeZdd<T>& n) const {
        OptimalSolution<T> sol(pi[num_jobs]);
        auto               m = std::min_element(n.list.begin(), n.list.end(),
                                  compare_sub_nodes<T>);
#ifndef NDEBUG
        auto weight = (*m)->weight;
#endif

        Label<SubNodeZdd<T>, T>* ptr_node = &((*m)->forward_label[0]);

        while (ptr_node->get_previous() != nullptr) {
            auto aux_prev_node = ptr_node->get_previous();
            auto aux_job = aux_prev_node->get_job();
            sol.C_max += aux_job->processing_time;
            sol.push_job_back(aux_job, aux_prev_node->get_weight(),
                              pi[aux_job->job]);
            ptr_node = aux_prev_node;
        }

        assert(sol.C_max == weight);

        return sol;
    }

    ForwardZddBase<T>(const ForwardZddBase<T>&) = default;
    ForwardZddBase<T>(ForwardZddBase<T>&&) noexcept = default;
    ForwardZddBase<T>& operator=(const ForwardZddBase<T>&) = default;
    ForwardZddBase<T>& operator=(ForwardZddBase<T>&&) noexcept = default;
};

template <typename T = double>
class ForwardZddCycle : public ForwardZddBase<T> {
    using ForwardZddBase<T>::pi;
    using ForwardZddBase<T>::num_jobs;

   public:
    ForwardZddCycle(T* _pi, int _num_jobs)
        : ForwardZddBase<T>(_pi, _num_jobs) {}

    explicit ForwardZddCycle(int _num_jobs) : ForwardZddBase<T>(_num_jobs) {}

    ForwardZddCycle() : ForwardZddBase<T>() {
        pi = nullptr;
        num_jobs = 0;
    }
    ~ForwardZddCycle<T>() = default;

    // ForwardZddCycle(const ForwardZddCycle<T>& src) {
    //     pi = src.pi;
    //     num_jobs = src.num_jobs;
    // }

    void initializenode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            if (it->weight == 0) {
                it->forward_label[0].update_solution(pi[num_jobs], nullptr,
                                                     false);
                it->forward_label[1].update_solution(
                    std::numeric_limits<double>::max() / 2, nullptr, false);
            } else {
                it->forward_label[0].update_solution(
                    std::numeric_limits<double>::max() / 2, nullptr, false);
                it->forward_label[1].update_solution(
                    std::numeric_limits<double>::max() / 2, nullptr, false);
            }
        }
    }

    void initializerootnode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->forward_label[0].get_f() = pi[num_jobs];
            it->forward_label[1].set_f(std::numeric_limits<double>::max() / 2);
        }
    }

    void evalNode(NodeZdd<T>& n) const override {
        Job* tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto& it : n.list) {
            int                            weight = it->weight;
            T                              g;
            std::shared_ptr<SubNodeZdd<T>> p0 = it->n;
            std::shared_ptr<SubNodeZdd<T>> p1 = it->y;
            double result = value_Fj(weight + tmp_j->processing_time, tmp_j) -
                            pi[tmp_j->job];

            /**
             * High edge calculation
             */
            Job* prev = it->forward_label[0].get_previous_job();
            Job* aux1 = p1->forward_label[0].get_previous_job();
            auto diff = (prev == nullptr)
                            ? true
                            : (value_diff_Fij(weight, tmp_j, prev) >= 0);

            if (prev != tmp_j && diff) {
                g = it->forward_label[0].get_f() + result;
                if (g < p1->forward_label[0].get_f()) {
                    if (aux1 != tmp_j) {
                        p1->forward_label[1].update_solution(
                            p1->forward_label[0]);
                    }
                    p1->forward_label[0].update_solution(
                        g, &(it->forward_label[0]), true);
                } else if ((g < p1->forward_label[1].get_f()) &&
                           (aux1 != tmp_j)) {
                    p1->forward_label[1].update_solution(
                        g, &(it->forward_label[0]), true);
                }
            } else {
                g = it->forward_label[1].get_f() + result;
                prev = it->forward_label[1].get_previous_job();
                diff = (prev == nullptr)
                           ? true
                           : (value_diff_Fij(weight, tmp_j, prev) >= 0);
                if (diff) {
                    if (g < p1->forward_label[0].get_f()) {
                        if (aux1 != tmp_j) {
                            p1->forward_label[1].update_solution(
                                p1->forward_label[0]);
                        }
                        p1->forward_label[0].update_solution(
                            g, &(it->forward_label[1]), true);
                    } else if ((g < p1->forward_label[1].get_f()) &&
                               (aux1 != tmp_j)) {
                        p1->forward_label[1].update_solution(
                            g, &(it->forward_label[1]), true);
                    }
                }
            }

            /**
             * Low edge calculation
             */
            aux1 = p0->forward_label[0].get_previous_job();
            if (it->forward_label[0].get_f() < p0->forward_label[0].get_f()) {
                if (prev != aux1) {
                    p0->forward_label[1].update_solution(p0->forward_label[0]);
                }
                p0->forward_label[0].update_solution(it->forward_label[0]);
                if (it->forward_label[1].get_f() <
                    p0->forward_label[1].get_f()) {
                    p0->forward_label[1].update_solution(it->forward_label[1]);
                }
            } else if ((it->forward_label[0].get_f() <
                        p0->forward_label[1].get_f()) &&
                       (aux1 != prev)) {
                p0->forward_label[1].update_solution(it->forward_label[0]);
            } else if ((it->forward_label[1].get_f() <
                        p0->forward_label[1].get_f())) {
                p0->forward_label[1].update_solution(it->forward_label[1]);
            }
        }
    }

    ForwardZddCycle<T>(const ForwardZddCycle<T>&) = default;
    ForwardZddCycle<T>(ForwardZddCycle<T>&&) noexcept = default;
    ForwardZddCycle<T>& operator=(const ForwardZddCycle<T>&) = default;
    ForwardZddCycle<T>& operator=(ForwardZddCycle<T>&&) noexcept = default;
};

template <typename T = double>
class ForwardZddSimple : public ForwardZddBase<T> {
    using ForwardZddBase<T>::pi;
    using ForwardZddBase<T>::num_jobs;

   public:
    ForwardZddSimple(T* _pi, int _num_jobs)
        : ForwardZddBase<T>(_pi, _num_jobs) {}

    explicit ForwardZddSimple(int _num_jobs) : ForwardZddBase<T>(_num_jobs) {}

    ForwardZddSimple() {
        pi = nullptr;
        num_jobs = 0;
    }

    ~ForwardZddSimple<T>() = default;

    // ForwardZddSimple(const ForwardZddSimple<T>& src) {
    //     pi = src.pi;
    //     num_jobs = src.num_jobs;
    // }

    void initializenode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            if (it->weight == 0) {
                it->forward_label[0].update_solution(pi[num_jobs], nullptr,
                                                     false);
            } else {
                it->forward_label[0].update_solution(
                    std::numeric_limits<double>::max() / 2, nullptr, false);
            }
        }
    }

    void initializerootnode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            // printf("test init %f\n", -pi[num_jobs]);
            it->forward_label[0].get_f() = pi[num_jobs];
        }
    }

    void initialize_pi(T* _pi) { pi = _pi; }

    void evalNode(NodeZdd<T>& n) const override {
        Job* tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto& it : n.list) {
            int                            weight = it->weight;
            T                              g;
            std::shared_ptr<SubNodeZdd<T>> p0 = it->n;
            std::shared_ptr<SubNodeZdd<T>> p1 = it->y;
            double result = value_Fj(weight + tmp_j->processing_time, tmp_j) -
                            pi[tmp_j->job];
            // printf("test result %f\n", result);

            /**
             * High edge calculation
             */
            g = it->forward_label[0].get_f() + result;
            if (g < p1->forward_label[0].get_f()) {
                p1->forward_label[0].update_solution(g, &(it->forward_label[0]),
                                                     true);
            }

            /**
             * Low edge calculation
             */
            if (it->forward_label[0].get_f() < p0->forward_label[0].get_f()) {
                p0->forward_label[0].update_solution(it->forward_label[0]);
            }
        }
    }
    ForwardZddSimple<T>(const ForwardZddSimple<T>&) = default;
    ForwardZddSimple<T>(ForwardZddSimple<T>&&) noexcept = default;
    ForwardZddSimple<T>& operator=(const ForwardZddSimple<T>&) = default;
    ForwardZddSimple<T>& operator=(ForwardZddSimple<T>&&) noexcept = default;
};

#endif  // DURATION_ZDD_HPP