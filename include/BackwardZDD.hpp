#ifndef BACKWARD_ZDD_HPP
#define BACKWARD_ZDD_HPP
#include <algorithm>
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"
#include "ZddNode.hpp"

template <typename T = double>
class BackwardZDDBase : public Eval<NodeZdd<T>, OptimalSolution<T>> {
    T*  pi{nullptr};
    int num_jobs{};

   protected:
   public:
    BackwardZDDBase(T* _pi, int _num_jobs) : pi(_pi), num_jobs(_num_jobs){};
    explicit BackwardZDDBase(int _num_jobs) : num_jobs(_num_jobs){};
    BackwardZDDBase() = default;
    BackwardZDDBase<T>(const BackwardZDDBase<T>&) = default;
    BackwardZDDBase<T>& operator=(const BackwardZDDBase<T>&) = default;
    BackwardZDDBase<T>(BackwardZDDBase<T>&&) = default;
    BackwardZDDBase<T>& operator=(BackwardZDDBase<T>&&) = default;
    ~BackwardZDDBase() = default;

    void initialize_pi(T* _pi) { pi = _pi; }

    virtual void initializenode(NodeZdd<T>& n) const = 0;
    virtual void initializerootnode(NodeZdd<T>& n) const = 0;
    virtual void evalNode(NodeZdd<T>& n) const = 0;
};

template <typename T = double>
class BackwardZddSimple : public BackwardZDDBase<T> {
   public:
    using BackwardZDDBase<T>::pi;
    using BackwardZDDBase<T>::num_jobs;

    BackwardZddSimple() : BackwardZDDBase<T>(){};
    BackwardZddSimple(T* _pi, int _num_jobs)
        : BackwardZDDBase<T>(_pi, _num_jobs){};
    explicit BackwardZddSimple(int _num_jobs) : BackwardZDDBase<T>(_num_jobs){};

    void evalNode(NodeZdd<T>& n) const override {
        Job* tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto& it : n.list) {
            auto weight = it->weight;
            auto p0 = it->n;
            auto p1 = it->y;
            auto result = value_Fj(weight + tmp_j->processing_time, tmp_j) -
                          pi[tmp_j->job];

            auto obj0 = p0->backward_label[0].get_f();
            auto obj1 = p1->backward_label[0].get_f() + result;

            if (obj0 > obj1) {
                it->backward_label[0].update_solution(obj1, nullptr, true);
            } else {
                it->backward_label[0].update_solution(obj0, nullptr, false);
            }
        }
    }

    void initializenode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
        }
    }

    void initializerootnode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].get_f() = pi[num_jobs];
        }
    }

    OptimalSolution<T> get_objective(NodeZdd<T>& n) const override {
        OptimalSolution<T> sol(pi[num_jobs]);

        auto m = std::min_element(n.list.begin(), n.list.end(),
                                  compare_sub_nodes<T>);
        std::shared_ptr<SubNodeZdd<T>> ptr_node = (*m);
        auto                           aux_job = n.get_job();

        // NodeZdd<T> *aux_node = &n;
        // Job *aux_job =  n.get_job();

        while (aux_job) {
            if (ptr_node->backward_label[0].get_high()) {
                sol.push_job_back(aux_job, pi[aux_job->job]);
                ptr_node = ptr_node->y;
                aux_job = ptr_node->get_job();
            } else {
                ptr_node = ptr_node->n;
                aux_job = ptr_node->get_job();
            }
        }

        return sol;
    }
};

template <typename T = double>
class BackwardZddCycle : public BackwardZDDBase<T> {
   public:
    using BackwardZDDBase<T>::pi;
    using BackwardZDDBase<T>::num_jobs;

    BackwardZddCycle() : BackwardZDDBase<T>(){};
    BackwardZddCycle(T* _pi, int _num_jobs)
        : BackwardZDDBase<T>(_pi, _num_jobs){};
    explicit BackwardZddCycle(int _num_jobs) : BackwardZDDBase<T>(_num_jobs){};

    void evalNode(NodeZdd<T>& n) const override {
        Job* tmp_j = n.get_job();

        for (auto& it : n.list) {
            int                           weight{it->get_weight()};
            std::shared_ptr<SubNodeZdd<>> p0{it->n};
            std::shared_ptr<SubNodeZdd<>> p1{it->y};
            T result{value_Fj(weight + tmp_j->processing_time, tmp_j) -
                     pi[tmp_j->job]};

            Job* prev_job{p1->backward_label[0].get_prev_job()};

            it->backward_label[0].update_label(&(p0->backward_label[0]));
            it->backward_label[1].update_label(&(p0->backward_label[1]));
            bool diff = bool_diff_Fij(weight, prev_job, tmp_j);
            bool diff1 = bool_diff_Fij(
                weight, p1->backward_label[0].get_prev_job(), tmp_j);

            if (prev_job != tmp_j && diff) {
                T obj1{p1->backward_label[0].get_f() + result};
                T obj2{p1->backward_label[1].get_f() + result};

                if (obj1 < it->backward_label[0].get_f()) {
                    if (tmp_j != it->backward_label[0].get_prev_job()) {
                        it->backward_label[1].update_label(
                            &(p0->backward_label[0]));
                    }

                    it->backward_label[0].update_label(&(p1->backward_label[0]),
                                                       obj1, true);
                } else if (obj1 < it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].get_prev_job() &&
                           diff1) {
                    it->backward_label[1].update_label(&(p1->backward_label[0]),
                                                       obj1, true);
                } else if (obj2 < it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].get_prev_job()) {
                    it->backward_label[1].update_label(&(p1->backward_label[1]),
                                                       obj2, true);
                }
            } else {
                T obj1 = p1->backward_label[1].get_f() + result;

                if (obj1 < it->backward_label[0].get_f()) {
                    if (tmp_j != it->backward_label[0].get_prev_job()) {
                        it->backward_label[1].update_label(
                            &(p0->backward_label[0]));
                    }

                    it->backward_label[0].update_label(&(p1->backward_label[1]),
                                                       obj1, true);
                } else if (obj1 < it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].get_prev_job()) {
                    it->backward_label[1].update_label(&(p1->backward_label[1]),
                                                       obj1, true);
                }
            }
        }
    }

    void initializenode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
        }
    }

    void initializerootnode(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].get_f() = pi[num_jobs];
        }
    }

    OptimalSolution<T> get_objective(NodeZdd<>& n) const override {
        OptimalSolution<T> sol(pi[num_jobs]);
        auto               m = std::min_element(n.list.begin(), n.list.end(),
                                  compare_sub_nodes<T>);
        Label<SubNodeZdd<T>, T>* aux_label = &((*m)->backward_label[0]);

        while (aux_label) {
            if (aux_label->get_high()) {
                Job* aux_job = aux_label->get_job();
                sol.push_job_back(aux_job, pi[aux_job->job]);
            }

            aux_label = aux_label->get_previous();
        }

        return sol;
    }
};

#endif  // BACKWARD_ZDD_HPP
