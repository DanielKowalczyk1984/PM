#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
#include <tdzdd/DdEval.hpp>
#include <OptimalSolution.hpp>
#include <node_duration.hpp>

template<typename E, typename T>
class BackwardBddBase : public tdzdd::DdEval<E, Node<T>, Optimal_Solution<T>> {
  protected:
    T *pi;
    int num_jobs;
  public:
    BackwardBddBase(T *_pi, int _num_jobs) : pi(_pi), num_jobs(_num_jobs) {};
    BackwardBddBase(int _num_jobs) : pi(nullptr), num_jobs(_num_jobs) {};
    BackwardBddBase(): pi(nullptr), num_jobs(0) {};
    BackwardBddBase(const BackwardBddBase<E, T> &src) : pi(src.pi),
        num_jobs(src.num_jobs) {};

    void initializepi(T *_pi) {
        pi = _pi;
    }

    virtual void initializenode(Node<T> &n) const  = 0;
    virtual void initializerootnode(Node<T> &n) const  = 0;
    virtual void evalNode(Node<T> &n) const = 0;

    Optimal_Solution<T> getValue(Node<T> const &n) {
        Optimal_Solution<T> sol;

        return sol;
    }
};

template<typename E, typename T>
class BackwardBddSimple : public BackwardBddBase<E, T> {
  public:
    using BackwardBddBase<E, T>::pi;
    using BackwardBddBase<E, T>::num_jobs;

    BackwardBddSimple() : BackwardBddBase<E, T>() {};
    BackwardBddSimple(T *_pi, int _num_jobs) : BackwardBddBase<E, T>(_pi,
                _num_jobs) {};
    explicit BackwardBddSimple(int _num_jobs) : BackwardBddBase<E, T>(_num_jobs) {};

    void evalNode(Node<T> &n) const override {
        Job *tmp_j = n.GetJob();
        int      weight = n.GetWeight();
        Node<T> *p0 = n.child[0];
        Node<T> *p1 = n.child[1];
        T result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

        T obj0 = p0->prev1.GetF();
        T obj1 = p1->prev1.GetF() + result;

        if (obj0 < obj1) {
            n.prev1.UpdateSolution(obj1, nullptr, true);
        } else {
            n.prev1.UpdateSolution(obj0, nullptr, false);
        }
    }

    void initializenode(Node<T> &n) const override {
        n.prev1.UpdateSolution(-DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(Node<T> &n) const override {
        n.prev1.f = pi[num_jobs];
    }

    Optimal_Solution<T> get_objective(Node<T> &n) const {
        Optimal_Solution<T> sol(pi[num_jobs]);

        Node<T> *aux_node = &n;
        Job *aux_job =  n.GetJob();

        while (aux_job) {
            if (aux_node->prev1.GetHigh()) {
                sol.push_job_back(aux_job, pi[aux_job->job]);
                aux_node = aux_node->child[1];
                aux_job = aux_node->GetJob();
            } else {
                aux_node = aux_node->child[0];
                aux_job = aux_node->GetJob();
            }

        }

        return sol;
    }
};

template<typename E, typename T>
class BackwardBddCycle : public BackwardBddBase<E, T> {
  public:
    using BackwardBddBase<E, T>::pi;
    using BackwardBddBase<E, T>::num_jobs;

    BackwardBddCycle() : BackwardBddBase<E, T>() {};
    BackwardBddCycle(T *_pi, int _num_jobs) : BackwardBddBase<E, T>(_pi,_num_jobs) {};
    explicit BackwardBddCycle(int _num_jobs) : BackwardBddBase<E, T>(_num_jobs) {};

    void evalNode(Node<T> &n) const override {
        Job *tmp_j = n.GetJob();
        int      weight = n.GetWeight();
        Node<T> *p0 = n.child[0];
        Node<T> *p1 = n.child[1];
        T result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

        Job *prev_job = p1->prev1.get_prev_job();

        if(prev_job != tmp_j) {
            T obj0 = p0->prev1.GetF();
            T obj1 = p1->prev1.GetF() + result;

            if(obj0 < obj1) {
                n.prev1.UpdateNode(obj1, tmp_j, true);

                prev_job = p0->prev1.get_prev_job();

                if(prev_job != tmp_j) {
                    n.prev2.UpdateNode(p0->prev1);
                } else {
                    n.prev2.UpdateNode(p0->prev2);
                }
            } else {
                n.prev1.UpdateNode(p0->prev1);
                obj0 = p0->prev2.GetF();

                if(obj0 >= obj1) {
                    n.prev2.UpdateNode(p0->prev2);
                } else {
                    if(tmp_j != n.prev1.get_prev_job()) {
                        n.prev2.UpdateNode(obj1, tmp_j, true);
                    } else {
                        obj1 = p1->prev2.GetF() + result;
                        if(obj0 >= obj1) {
                            n.prev2.UpdateNode(p0->prev2);
                        } else {
                            n.prev2.UpdateNode(obj1, tmp_j, true);
                        }
                    }
                }
            }
        } else {
            T obj0 = p0->prev1.GetF();
            T obj1 = p1->prev2.GetF() + result;

            if(obj0 < obj1) {
                n.prev1.UpdateNode(obj1, tmp_j, true);
                prev_job = p0->prev1.get_prev_job();
                if(prev_job != tmp_j) {
                    n.prev2.UpdateNode(p0->prev1);
                } else {
                    n.prev2.UpdateNode(p0->prev2);
                }
            } else {
                n.prev1.UpdateNode(p0->prev1);
                n.prev2.UpdateNode(p0->prev2);
            }
        }
    }

    void initializenode(Node<T> &n) const override {
        n.prev1.UpdateSolution(-DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(Node<T> &n) const override {
        n.prev1.f = pi[num_jobs];
    }

    Optimal_Solution<T> get_objective(Node<T> &n) const {
        Optimal_Solution<T> sol(pi[num_jobs]);

        // Node<T> *aux_node = &n;
        // Job *aux_job =  n.GetJob();

        // while (aux_job) {
        //     if (aux_node->prev1.GetHigh()) {
        //         sol.push_job_back(aux_job, pi[aux_job->job]);
        //         aux_node = aux_node->child[1];
        //         aux_job = aux_node->GetJob();
        //     } else {
        //         aux_node = aux_node->child[0];
        //         aux_job = aux_node->GetJob();
        //     }

        // }

        return sol;
    }
};

#endif // BACKWARD_BDD_HPP
