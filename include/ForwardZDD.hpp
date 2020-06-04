#ifndef DURATION_ZDD_HPP
#define DURATION_ZDD_HPP
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"
#include "ZddNode.hpp"
#include <vector>
#include <algorithm>
using namespace std;


template<typename E, typename T> class ForwardZddBase : public 
    Eval<E, NodeZdd<T>, OptimalSolution<T>> {
protected:
    T *pi;
    int num_jobs;
public:
    ForwardZddBase(T *_pi, int _num_jobs): pi(_pi), num_jobs(_num_jobs) {
    }

    ForwardZddBase(int _num_jobs)
    :  pi(nullptr), num_jobs(_num_jobs){
    }

    ForwardZddBase(){
        pi = nullptr;
        num_jobs = 0;
    }

    ForwardZddBase(const ForwardZddBase<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initialize_pi(T *_pi){
        pi = _pi;
    }

    virtual void initializenode(NodeZdd<T>& n) const  = 0;

    virtual void initializerootnode(NodeZdd<T>& n) const  = 0;

    virtual void evalNode(NodeZdd<T>& n) const = 0;

    OptimalSolution<T> get_objective(NodeZdd<T> &n) const {
        OptimalSolution<T> sol(pi[num_jobs]);
        auto m =  std::min_element(n.list.begin(), n.list.end(),compare_sub_nodes<T>);
        #ifndef NDEBUG
        auto weight = (*m)->weight;
        #endif

        Label<SubNodeZdd<T>,T> *ptr_node = &((*m)->forward_label[0]);

        while(ptr_node->get_previous() != nullptr) {
            auto aux_prev_node = ptr_node->get_previous();
            auto aux_job = aux_prev_node->get_job();
            sol.C_max += aux_job->processing_time;
            sol.push_job_back(aux_job, aux_prev_node->get_weight(), pi[aux_job->job]);
            ptr_node = aux_prev_node;
        }

        assert(sol.C_max == weight);

        return sol;
    }

};

template<typename E, typename T> class ForwardZddCycle : public ForwardZddBase<E, T> {

  public:
    using ForwardZddBase<E, T>::pi;
    using ForwardZddBase<E, T>::num_jobs;
    ForwardZddCycle(T *_pi, int _num_jobs): ForwardZddBase<E, T>(_pi, _num_jobs) {
    }

    explicit ForwardZddCycle(int _num_jobs) : ForwardZddBase<E, T>(_num_jobs){
    }

    ForwardZddCycle(): ForwardZddBase<E, T>(){
        pi = nullptr;
        num_jobs = 0;
    }

    ForwardZddCycle(const ForwardZddCycle<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initializenode(NodeZdd<T>& n) const override {
        for (auto &it: n.list) {
            if(it->weight == 0) {
                it->forward_label[0].update_solution(pi[num_jobs], nullptr, false);
                it->forward_label[1].update_solution(DBL_MAX/2, nullptr, false);
            } else {
                it->forward_label[0].update_solution(DBL_MAX/2, nullptr, false);
                it->forward_label[1].update_solution(DBL_MAX/2, nullptr, false);
            }
        }
    }

    void initializerootnode(NodeZdd<T> &n) const override {
        for(auto &it: n.list){
            it->forward_label[0].f = pi[num_jobs];
            it->forward_label[1].set_f(DBL_MAX/2);
        }
    }

    void evalNode(NodeZdd<T> &n) const override
    {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto &it : n.list) {
            int      weight = it->weight;
            T g;
            std::shared_ptr<SubNodeZdd<T>> p0 = it->n;
            std::shared_ptr<SubNodeZdd<T>> p1 = it->y;
            double result = value_Fj(weight + tmp_j->processing_time, tmp_j) - pi[tmp_j->job];

            /**
             * High edge calculation
             */
            Job *prev = it->forward_label[0].get_previous_job();
            Job *aux1 = p1->forward_label[0].get_previous_job();
            double diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );

            if(prev != tmp_j && diff) {
                g = it->forward_label[0].get_f() + result;
                if(g < p1->forward_label[0].get_f()) {
                    if(aux1 != tmp_j) {
                        p1->forward_label[1].update_solution(p1->forward_label[0]);
                    }
                    p1->forward_label[0].update_solution(g, &(it->forward_label[0]), true);
                } else if ((g < p1->forward_label[1].get_f()) && (aux1 != tmp_j)) {
                    p1->forward_label[1].update_solution(g, &(it->forward_label[0]), true);
                }
            } else  {
                g = it->forward_label[1].get_f() + result;
                prev = it->forward_label[1].get_previous_job();
                diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );
                if(diff) {
                    if(g < p1->forward_label[0].get_f()) {
                        if(aux1 != tmp_j) {
                            p1->forward_label[1].update_solution(p1->forward_label[0]);
                        }
                        p1->forward_label[0].update_solution(g, &(it->forward_label[1]), true);
                    } else if ((g < p1->forward_label[1].get_f()) && (aux1 != tmp_j)) {
                        p1->forward_label[1].update_solution(g, &(it->forward_label[1]), true);
                    }
                }
            }

            /**
             * Low edge calculation
             */
            aux1 = p0->forward_label[0].get_previous_job();
            if(it->forward_label[0].get_f() < p0->forward_label[0].get_f()) {
                if(prev != aux1) {
                    p0->forward_label[1].update_solution(p0->forward_label[0]);
                }
                p0->forward_label[0].update_solution(it->forward_label[0]);
                if(it->forward_label[1].get_f() < p0->forward_label[1].get_f()) {
                    p0->forward_label[1].update_solution(it->forward_label[1]);
                }
            } else if ((it->forward_label[0].get_f() < p0->forward_label[1].get_f()) && (aux1 != prev)){
                p0->forward_label[1].update_solution(it->forward_label[0]);
            } else if ((it->forward_label[1].get_f() < p0->forward_label[1].get_f())) {
                p0->forward_label[1].update_solution(it->forward_label[1]);
            }
        }
    }

    OptimalSolution<T> getValue(NodeZdd<T> const &n) {
        OptimalSolution<T> sol;

        return sol;
    }
};

template<typename E, typename T> class ForwardZddSimple : public ForwardZddBase<E, T> {
  public:
    using ForwardZddBase<E, T>::pi;
    using ForwardZddBase<E, T>::num_jobs;
    ForwardZddSimple(T *_pi, int _num_jobs): ForwardZddBase<E, T>(_pi, _num_jobs) {
    }

    explicit ForwardZddSimple(int _num_jobs)
    :  ForwardZddBase<E, T> (_num_jobs){
    }

    ForwardZddSimple(){
        pi = nullptr;
        num_jobs = 0;
    }

    ForwardZddSimple(const ForwardZddSimple<E, T> &src) {
        pi = src.pi;
        num_jobs = src.num_jobs;
    }

    void initializenode(NodeZdd<T>& n) const override {
        for (auto &it: n.list) {
            if(it->weight == 0) {
                it->forward_label[0].update_solution(pi[num_jobs], nullptr, false);
            } else {
                it->forward_label[0].update_solution(DBL_MAX/2, nullptr, false);
            }
        }
    }

    void initializerootnode(NodeZdd<T> &n) const override {
        for(auto &it: n.list){
            // printf("test init %f\n", -pi[num_jobs]);
            it->forward_label[0].f = pi[num_jobs];
        }
    }

    void initialize_pi(T *_pi){
        pi = _pi;
    }

    void evalNode(NodeZdd<T> &n) const override {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto &it : n.list) {
            int      weight = it->weight;
            T g;
            std::shared_ptr<SubNodeZdd<T>> p0 = it->n;
            std::shared_ptr<SubNodeZdd<T>> p1 = it->y;
            double result = value_Fj(weight + tmp_j->processing_time, tmp_j) - pi[tmp_j->job];
            // printf("test result %f\n", result);

            /**
             * High edge calculation
             */
            g = it->forward_label[0].get_f() + result;
            if(g < p1->forward_label[0].get_f()) {
                p1->forward_label[0].update_solution(g, &(it->forward_label[0]), true);
            }

            /**
             * Low edge calculation
             */
            if(it->forward_label[0].get_f() < p0->forward_label[0].get_f()) {
                p0->forward_label[0].update_solution(it->forward_label[0]);
            }
        }
    }

    OptimalSolution<T> getValue(NodeZdd<T> const &n) {
        OptimalSolution<T> sol;
        return sol;
    }
};

#endif // DURATION_ZDD_HPP