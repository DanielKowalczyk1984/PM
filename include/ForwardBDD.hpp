#ifndef FORWARD_BDD_HPP
#define FORWARD_BDD_HPP
// #include <tdzdd/DdEval.hpp>
#include "NodeBddEval.hpp"
#include <OptimalSolution.hpp>
#include <node_duration.hpp>


template<typename E, typename T> class ForwardBddBase : public 
    NodeBddEval<E, Node<T>, OptimalSolution<T>> {
protected:
    T *pi;
    int num_jobs;
public:
    ForwardBddBase(T *_pi, int _num_jobs): pi(_pi), num_jobs(_num_jobs) {
    }

    ForwardBddBase(int _num_jobs)
    :  pi(nullptr), num_jobs(_num_jobs){
    }

    ForwardBddBase(){
        pi = nullptr;
        num_jobs = 0;
    }

    ForwardBddBase(const ForwardBddBase<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initializepi(T *_pi){
        pi = _pi;
    }

    virtual void initializenode(Node<T>& n) const = 0;

    virtual void initializerootnode(Node<T>& n) const  = 0;

    virtual void evalNode(Node<T>& n) const = 0;

    OptimalSolution<T> get_objective(Node<T> &n) const {
        OptimalSolution<T> sol(-pi[num_jobs]);
        Label<T> *ptr_node = &(n.forward_label1);

        while(ptr_node->get_previous() != nullptr) {
            Label<T> *aux_prev_node = ptr_node->get_previous();
            Job *aux_job = aux_prev_node->get_job();
            sol.C_max += aux_job->processing_time;
            sol.push_job_back(aux_job, aux_prev_node->get_weight(), pi[aux_job->job]);
            ptr_node = aux_prev_node;
        }

        return sol;
    }

    OptimalSolution<T> getValue(Node<T> const &n){
        OptimalSolution<T> sol;

        return sol;
    }
};


template<typename E, typename T> class ForwardBddCycle : public ForwardBddBase<E, T> {
  public:
    using ForwardBddBase<E, T>::pi;
    using ForwardBddBase<E, T>::num_jobs;

    ForwardBddCycle(T *_pi, int _num_jobs) : ForwardBddBase<E, T>(_pi, _num_jobs) {}

    explicit ForwardBddCycle(int _num_jobs) : ForwardBddBase<E, T>(_num_jobs) {}

    ForwardBddCycle(): ForwardBddBase<E, T>(){
        pi = nullptr;
        num_jobs = 0;
    }

    ForwardBddCycle(const ForwardBddCycle<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initializenode(Node<T>& n) const override {
        if(n.get_weight() == 0) {
            n.forward_label1.update_solution(-pi[num_jobs], nullptr, false);
            n.forward_label2.update_solution(-DBL_MAX/2, nullptr, false);
        } else {
            n.forward_label1.update_solution(-DBL_MAX/2, nullptr, false);
            n.forward_label2.update_solution(-DBL_MAX/2, nullptr, false);
        }
    }

    void initializerootnode(Node<T> &n) const override {
        n.forward_label1.f = -pi[num_jobs];
        n.forward_label2.set_f(-DBL_MAX/2);
    }

    void evalNode(Node<T> &n) const override
    {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);
        double result;
        bool diff;

        int      weight = n.get_weight();
        T g;
        Node<T>* p0 = n.child[0];
        Node<T>* p1 = n.child[1];
        result = - value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job];

        /**
         * High edge calculation
         */
        Job *prev = n.forward_label1.get_previous_job();
        Job *aux1 = p1->forward_label1.get_previous_job();
        diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );

        if(prev != tmp_j && diff) {
            g = n.forward_label1.get_f() + result;
            if(g > p1->forward_label1.get_f()) {
                if(aux1 != tmp_j) {
                    p1->forward_label2.update_solution(p1->forward_label1);
                }
                p1->forward_label1.update_solution(g, &(n.forward_label1), true);
            } else if ((g > p1->forward_label2.get_f()) && (aux1 != tmp_j)) {
                p1->forward_label2.update_solution(g, &(n.forward_label1), true);
            }
        } else  {
            g = n.forward_label2.get_f() + result;
            prev = n.forward_label2.get_previous_job();
            diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );

            if(diff) {
                if(g > p1->forward_label1.get_f()) {
                    if(aux1 != tmp_j) {
                        p1->forward_label2.update_solution(p1->forward_label1);
                    }
                    p1->forward_label1.update_solution(g, &(n.forward_label2), true);
                } else if ((g > p1->forward_label2.get_f()) && (aux1 != tmp_j)) {
                    p1->forward_label2.update_solution(g, &(n.forward_label2), true);
                }
            }
        }

        /**
         * Low edge calculation
         */
        aux1 = p0->forward_label1.get_previous_job();
        if(n.forward_label1.get_f() > p0->forward_label1.get_f()) {
            if(prev != aux1) {
                p0->forward_label2.update_solution(p0->forward_label1);
            }
            p0->forward_label1.update_solution(n.forward_label1);
            if(n.forward_label2.get_f() > p0->forward_label2.get_f()) {
                p0->forward_label2.update_solution(n.forward_label2);
            }
        } else if ((n.forward_label1.get_f() > p0->forward_label2.get_f()) && (aux1 != prev)){
            p0->forward_label2.update_solution(n.forward_label1);
        } else if ((n.forward_label2.get_f() > p0->forward_label2.get_f())) {
            p0->forward_label2.update_solution(n.forward_label2);
        }
    }
};


template<typename E, typename T> class ForwardBddSimple : public ForwardBddBase<E, T> {
  public:
    using ForwardBddBase<E, T>::pi;
    using ForwardBddBase<E, T>::num_jobs;
    ForwardBddSimple(T *_pi, int _num_jobs): ForwardBddBase<E, T>(_pi, _num_jobs) {
    }

    explicit ForwardBddSimple(int _num_jobs)
    :  ForwardBddBase<E, T>(_num_jobs){
    }

    ForwardBddSimple(){
        pi = nullptr;
        num_jobs = 0;
    }

    ForwardBddSimple(const ForwardBddSimple<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initializenode(Node<T>& n) const override {
        if(n.get_weight() == 0) {
            n.forward_label1.update_solution(-pi[num_jobs], nullptr, false);
        } else {
            n.forward_label1.update_solution(-DBL_MAX/2, nullptr, false);
        }
    }

    void initializerootnode(Node<T> &n) const override {
        n.forward_label1.f = -pi[num_jobs];
    }

    void initializepi(T *_pi){
        pi = _pi;
    }

    void evalNode(Node<T> &n) const override {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);
        T result;

        int      weight = n.get_weight();
        T g;
        Node<T>* p0 = n.child[0];
        Node<T>* p1 = n.child[1];
        result = - value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job];

        /**
         * High edge calculation
         */
        g = n.forward_label1.get_f() + result;
        if(g > p1->forward_label1.get_f()) {
            p1->forward_label1.update_solution(g, &(n.forward_label1), true);
        }

        /**
         * Low edge calculation
         */
        if(n.forward_label1.get_f() > p0->forward_label1.get_f()) {
            p0->forward_label1.update_solution(n.forward_label1);
        }
    }
};

#endif // FORWARD_BDD_HPP
