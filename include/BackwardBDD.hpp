#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
// #include <tdzdd/DdEval.hpp>
#include "NodeBddEval.hpp"
#include <OptimalSolution.hpp>
#include <node_duration.hpp>

template<typename E, typename T>
class BackwardBddBase : public NodeBddEval<E, Node<T>, OptimalSolution<T> > {
  protected:
    T *pi;
    int num_jobs;
  public:
    BackwardBddBase(T *_pi, int _num_jobs) : pi(_pi), num_jobs(_num_jobs) {
    };

    BackwardBddBase(int _num_jobs) : pi(nullptr), num_jobs(_num_jobs) {

    };

    BackwardBddBase() : pi(nullptr), num_jobs(0) {
    };

    BackwardBddBase(const BackwardBddBase<E, T> &src) : pi(src.pi),
        num_jobs(src.num_jobs) {
    };

    void initializepi(T *_pi) {
        pi = _pi;
    }

    virtual void initializenode(Node<T> &n) const  = 0;
    virtual void initializerootnode(Node<T> &n) const  = 0;
    virtual void evalNode(Node<T> &n) const = 0;
    virtual OptimalSolution<T> getValue(Node<T> const &n) = 0;
};

template<typename E, typename T>
class BackwardBddSimple : public BackwardBddBase<E, T> {
  public:
    using BackwardBddBase<E, T>::pi;
    using BackwardBddBase<E, T>::num_jobs;

    BackwardBddSimple() : BackwardBddBase<E, T>() {
    };
    BackwardBddSimple(T *_pi, int _num_jobs) : BackwardBddBase<E, T>(_pi,
                _num_jobs) {
    };
    explicit BackwardBddSimple(int _num_jobs) : BackwardBddBase<E, T>(_num_jobs) {
    };

    void evalNode(Node<T> &n) const override {
        Job *tmp_j = n.get_job();
        int weight = n.get_weight();
        Node<T> *p0 = n.child[0];
        Node<T> *p1 = n.child[1];
        T result = -value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job];

        T obj0 = p0->backward_label1.get_f();
        T obj1 = p1->backward_label1.get_f() + result;

        if (obj0 < obj1) {
            n.backward_label1.update_solution(obj1, nullptr, true);
        } else {
            n.backward_label1.update_solution(obj0, nullptr, false);
        }

    }

    void initializenode(Node<T> &n) const override {
        n.backward_label1.update_solution(-DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(Node<T> &n) const override {
        n.backward_label1.f = -pi[num_jobs];
    }

    OptimalSolution<T> get_objective(Node<T> &n) const {
        OptimalSolution<T> sol(-pi[num_jobs]);

        Node<T> *aux_node = &n;
        Job *aux_job =  n.get_job();

        while (aux_job) {
            if (aux_node->backward_label1.get_high()) {
                sol.push_job_back(aux_job, pi[aux_job->job]);
                aux_node = aux_node->child[1];
                aux_job = aux_node->get_job();
            } else {
                aux_node = aux_node->child[0];
                aux_job = aux_node->get_job();
            }

        }

        return sol;
    }

    OptimalSolution<T> getValue(Node<T> const &n) override {
        OptimalSolution<T> sol;
        return sol;
    }

};



template<typename E, typename T>
class BackwardBddCycle : public BackwardBddBase<E, T> {
  public:
    using BackwardBddBase<E, T>::pi;
    using BackwardBddBase<E, T>::num_jobs;

    BackwardBddCycle() : BackwardBddBase<E, T>() {
    };
    BackwardBddCycle(T *_pi, int _num_jobs) : BackwardBddBase<E, T>(_pi,
                _num_jobs) {
    };
    explicit BackwardBddCycle(int _num_jobs) : BackwardBddBase<E, T>(_num_jobs) {
    };

    void evalNode(Node<T> &n) const override {
        Job *tmp_j = n.get_job();
        int weight{n.get_weight()};
        Node<T> *p0  {n.child[0]};
        Node<T> *p1  {n.child[1]};
        T result { -value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job]};

        Job *prev_job{p1->backward_label1.get_prev_job()};

        n.backward_label1.update_label(&(p0->backward_label1));
        n.backward_label2.update_label(&(p0->backward_label2));
        bool diff =  bool_diff_Fij(weight, prev_job, tmp_j);
        bool diff1 = bool_diff_Fij(weight, p1->backward_label1.get_prev_job(), tmp_j);
        bool diff2 = bool_diff_Fij(weight, p1->backward_label2.get_prev_job(), tmp_j);

        if (prev_job != tmp_j && diff) {
            T obj1 {p1->backward_label1.get_f() + result};
            T obj2 {p1->backward_label2.get_f() + result};

            if (obj1 > n.backward_label1.get_f()) {
                if (tmp_j != n.backward_label1.get_prev_job()) {
                    n.backward_label2.update_label(&(p0->backward_label1));
                }

                n.backward_label1.update_label(&(p1->backward_label1), obj1, true);
            } else if (obj1 > n.backward_label2.get_f() &&
                       tmp_j != n.backward_label1.get_prev_job() && diff1) {
                n.backward_label2.update_label(&(p1->backward_label1), obj1, true);
            } else if (obj2 > n.backward_label2.get_f() &&
                       tmp_j != n.backward_label1.get_prev_job() && diff2) {
                n.backward_label2.update_label(&(p1->backward_label2), obj2, true);
            }
        } else {
            T obj1 = p1->backward_label2.get_f() + result;

            if (obj1 > n.backward_label1.get_f() && diff2) {
                if (tmp_j != n.backward_label1.get_prev_job()) {
                    n.backward_label2.update_label(&(p0->backward_label1));
                }

                n.backward_label1.update_label(&(p1->backward_label2), obj1, true);
            } else if (obj1 > n.backward_label2.get_f() &&
                       tmp_j != n.backward_label1.get_prev_job() && diff2) {
                n.backward_label2.update_label(&(p1->backward_label2), obj1, true);
            }
        }
    }

    void initializenode(Node<T> &n) const override {
        n.backward_label1.update_solution(-DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(Node<T> &n) const override {
        n.backward_label1.f = -pi[num_jobs];
    }

    OptimalSolution<T> get_objective(Node<T> &n) const {
        OptimalSolution<T> sol(-pi[num_jobs]);
        Label<T> *aux_label = &(n.backward_label1);

        while (aux_label) {
            if (aux_label->get_high()) {
                Job *aux_job = aux_label->get_job();
                sol.push_job_back(aux_job, pi[aux_job->job]);
            }

            aux_label = aux_label->get_previous();
        }

        return sol;
    }

    OptimalSolution<T> getValue(Node<T> const &n) override {
        OptimalSolution<T> sol;
        return sol;
    }
};

#endif // BACKWARD_BDD_HPP
