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

        T obj0 = p0->state1.GetF();
        T obj1 = p1->state1.GetF() + result;

        if (obj0 < obj1) {
            n.state1.UpdateSolution(obj1, nullptr, true);
        } else {
            n.state1.UpdateSolution(obj0, nullptr, false);
        }
    }

    void initializenode(Node<T> &n) const override {
        n.state1.UpdateSolution(-DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(Node<T> &n) const override {
        n.state1.f = -pi[num_jobs];
    }

    Optimal_Solution<T> get_objective(Node<T> &n) const {
        Optimal_Solution<T> sol(-pi[num_jobs]);

        Node<T> *aux_node = &n;
        Job *aux_job =  n.GetJob();

        while (aux_job) {
            if (aux_node->state1.GetHigh()) {
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
        int      weight{n.GetWeight()};
        Node<T> *p0  {n.child[0]};
        Node<T> *p1  {n.child[1]};
        T result {- value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job]};

        Job *prev_job{p1->state1.get_prev_job()};

        n.state1.UpdateNode(&(p0->state1));
        n.state2.UpdateNode(&(p0->state2));
        bool diff =  bool_diff_Fij(weight, prev_job, tmp_j);
        bool diff1 = bool_diff_Fij(weight, p1->state1.get_prev_job(), tmp_j); 
        bool diff2 = bool_diff_Fij(weight, p1->state2.get_prev_job(), tmp_j); 

        if(prev_job != tmp_j && diff) {
            T obj1 {p1->state1.GetF() + result};
            T obj2 {p1->state2.GetF() + result};

            if(obj1 > n.state1.GetF()) {
                if(tmp_j != n.state1.get_prev_job()) {
                    n.state2.UpdateNode(&(p0->state1));
                }
                n.state1.UpdateNode(&(p1->state1), obj1, true);
            } else if (obj1 > n.state2.GetF() && tmp_j != n.state1.get_prev_job() && diff1){
                n.state2.UpdateNode(&(p1->state1), obj1, true);
            } else if (obj2 > n.state2.GetF() && tmp_j != n.state1.get_prev_job() && diff2) {
                n.state2.UpdateNode(&(p1->state2), obj2, true);
            }
        } else {
            T obj1 = p1->state2.GetF() + result;
            if(obj1 > n.state1.GetF()) {
                if(tmp_j != n.state1.get_prev_job()) {
                    n.state2.UpdateNode(&(p0->state1));
                }
                n.state1.UpdateNode(&(p1->state2), obj1, true);
            } else if (obj1 > n.state2.GetF() && tmp_j != n.state1.get_prev_job()){
                n.state2.UpdateNode(&(p1->state2), obj1, true);
            }
        }
    }

    void initializenode(Node<T> &n) const override {
        n.state1.UpdateSolution(-DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(Node<T> &n) const override {
        n.state1.f = -pi[num_jobs];
    }

    Optimal_Solution<T> get_objective(Node<T> &n) const {
        Optimal_Solution<T> sol(-pi[num_jobs]);

        Label<T> *aux_node = &n.state1;
        Job *aux_job =  aux_node->GetJob();

        while (aux_job) {
            if (aux_node->GetHigh()) {
                sol.push_job_back(aux_job, pi[aux_job->job]);
                aux_node = aux_node->GetPrev();
                aux_job = aux_node->GetJob();
            } else {
                aux_node = aux_node->GetPrev();
                aux_job = aux_node->GetJob();
            }

        }

        return sol;
    }
};

#endif // BACKWARD_BDD_HPP
