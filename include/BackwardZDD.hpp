#ifndef BACKWARD_ZDD_HPP
#define BACKWARD_ZDD_HPP
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"
#include "ZddNode.hpp"
#include <algorithm>

template <typename E, typename T>
class BackwardZDDBase : public Eval<E, NodeZdd<T>, OptimalSolution<T>> {
 protected:
  T* pi;
  int num_jobs;

 public:
  BackwardZDDBase(T* _pi, int _num_jobs) : pi(_pi), num_jobs(_num_jobs){};
  explicit BackwardZDDBase(int _num_jobs) : pi(nullptr), num_jobs(_num_jobs){};
  BackwardZDDBase() : pi(nullptr), num_jobs(0){};
  ~BackwardZDDBase() {};

  void initialize_pi(T *_pi) {
      pi = _pi;
  }

  virtual void initializenode(NodeZdd<T> &n) const  = 0;
  virtual void initializerootnode(NodeZdd<T> &n) const  = 0;
  virtual void evalNode(NodeZdd<T> &n) const = 0;
  virtual OptimalSolution<T> getValue(NodeZdd<T> const &n) = 0;

};


template<typename E, typename T>
class BackwardZddSimple : public BackwardZDDBase<E, T> {
  public:
    using BackwardZDDBase<E, T>::pi;
    using BackwardZDDBase<E, T>::num_jobs;

    BackwardZddSimple() : BackwardZDDBase<E, T>() {
    };
    BackwardZddSimple(T *_pi, int _num_jobs) : BackwardZDDBase<E, T>(_pi,
                _num_jobs) {
    };
    explicit BackwardZddSimple(int _num_jobs) : BackwardZDDBase<E, T>(_num_jobs) {
    };

    void evalNode(NodeZdd<T> &n) const override {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for(auto &it : n.list) {
            auto weight = it->weight;
            auto p0 = it->n;
            auto p1 = it->y;
            auto result = -value_Fj(weight + tmp_j->processing_time,tmp_j) + pi[tmp_j->job];

            auto obj0 = p0->backward_label[0].get_f();
            auto obj1 = p1->backward_label[0].get_f() + result;

            if(obj0 < obj1) {
                it->backward_label[0].update_solution(obj1,nullptr,true);
            } else {
                it->backward_label[0].update_solution(obj0,nullptr,false);
            }

        }
    }

    void initializenode(NodeZdd<T> &n) const override {
        for(auto &it : n.list) {
            it->backward_label[0].update_solution(-DBL_MAX / 2, nullptr, false);
        }
    }

    void initializerootnode(NodeZdd<T> &n) const override {
        for(auto &it : n.list) {
            it->backward_label[0].f = -pi[num_jobs];
        }
    }

    OptimalSolution<T> get_objective(NodeZdd<T> &n) const {
        OptimalSolution<T> sol(-pi[num_jobs]);

        auto m =  std::max_element(n.list.begin(),n.list.end(),compare_sub_nodes<T>);
        std::shared_ptr<SubNodeZdd<T>> ptr_node = (*m);
        auto aux_job = n.get_job();

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

    OptimalSolution<T> getValue(NodeZdd<T> const &n) override {
        OptimalSolution<T> sol;
        return sol;
    }

};

template<typename E, typename T>
class BackwardZddCycle : public BackwardZDDBase<E, T> {
  public:
    using BackwardZDDBase<E, T>::pi;
    using BackwardZDDBase<E, T>::num_jobs;

    BackwardZddCycle() : BackwardZDDBase<E, T>() {
    };
    BackwardZddCycle(T *_pi, int _num_jobs) : BackwardZDDBase<E, T>(_pi,
                _num_jobs) {
    };
    explicit BackwardZddCycle(int _num_jobs) : BackwardZDDBase<E, T>(_num_jobs) {
    };

    void evalNode(NodeZdd<T> &n) const override {
        Job *tmp_j = n.get_job();

        for(auto &it : n.list) {
            int weight{it->get_weight()};
            std::shared_ptr<SubNodeZdd<>> p0  {it->n};
            std::shared_ptr<SubNodeZdd<>> p1  {it->y};
            T result { -value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job]};

            Job *prev_job{p1->backward_label[0].get_prev_job()};

            it->backward_label[0].update_label(&(p0->backward_label[0]));
            it->backward_label[1].update_label(&(p0->backward_label[1]));
            bool diff =  bool_diff_Fij(weight, prev_job, tmp_j);
            bool diff1 = bool_diff_Fij(weight, p1->backward_label[0].get_prev_job(), tmp_j);
            bool diff2 = bool_diff_Fij(weight, p1->backward_label[1].get_prev_job(), tmp_j);

            if (prev_job != tmp_j && diff) {
                T obj1 {p1->backward_label[0].get_f() + result};
                T obj2 {p1->backward_label[1].get_f() + result};

                if (obj1 > it->backward_label[0].get_f()) {
                    if (tmp_j != it->backward_label[0].get_prev_job()) {
                        it->backward_label[1].update_label(&(p0->backward_label[0]));
                    }

                    it->backward_label[0].update_label(&(p1->backward_label[0]), obj1, true);
                } else if (obj1 > it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].get_prev_job() && diff1) {
                    it->backward_label[1].update_label(&(p1->backward_label[0]), obj1, true);
                } else if (obj2 > it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].get_prev_job() && diff2) {
                    it->backward_label[1].update_label(&(p1->backward_label[1]), obj2, true);
                }
            } else {
                T obj1 = p1->backward_label[1].get_f() + result;

                if (obj1 > it->backward_label[0].get_f() && diff2) {
                    if (tmp_j != it->backward_label[0].get_prev_job()) {
                        it->backward_label[1].update_label(&(p0->backward_label[0]));
                    }

                    it->backward_label[0].update_label(&(p1->backward_label[1]), obj1, true);
                } else if (obj1 > it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].get_prev_job() && diff2) {
                    it->backward_label[1].update_label(&(p1->backward_label[1]), obj1, true);
                }
            }

        }

    }

    void initializenode(NodeZdd<T> &n) const override {
        for(auto &it:n.list) {
            it->backward_label[0].update_solution(-DBL_MAX / 2, nullptr, false);
        }
    }

    void initializerootnode(NodeZdd<T> &n) const override {
        for(auto &it : n.list) {
            it->backward_label[0].f = -pi[num_jobs];
        }
    }

    OptimalSolution<T> get_objective(NodeZdd<> &n) const override {
        OptimalSolution<T> sol(-pi[num_jobs]);
        auto m = std::max_element(n.list.begin(), n.list.end(), compare_sub_nodes<T>);
        Label<SubNodeZdd<T>,T> *aux_label = &((*m)->backward_label[0]);

        while (aux_label) {
            if (aux_label->get_high()) {
                Job *aux_job = aux_label->get_job();
                sol.push_job_back(aux_job, pi[aux_job->job]);
            }

            aux_label = aux_label->get_previous();
        }

        return sol;
    }

    OptimalSolution<T> getValue(NodeZdd<T> const &n) override {
        OptimalSolution<T> sol;
        return sol;
    }
};

#endif  // BACKWARD_ZDD_HPP
