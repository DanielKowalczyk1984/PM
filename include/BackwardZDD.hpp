#ifndef BACKWARD_ZDD_HPP
#define BACKWARD_ZDD_HPP
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"
#include "ZddNode.hpp"
#include <algorithm>


// template<typename T>
// bool my_compare(const std::shared_ptr<SubNodeZdd<T>> &lhs, const std::shared_ptr<SubNodeZdd<T>> &rhs){
//     return *lhs < *rhs;
// }

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

  void initializepi(T *_pi) {
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

        auto m =  std::max_element(n.list.begin(),n.list.end(),my_compare<T>);
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

#endif  // BACKWARD_ZDD_HPP
