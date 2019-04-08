#ifndef DURATION_ZDD_HPP
#define DURATION_ZDD_HPP
#include <tdzdd/DdEval.hpp>
#include <OptimalSolution.hpp>
#include <node_duration.hpp>
#include <vector>
#include <algorithm>
using namespace std;

template<typename T>
class ForwardZddNode {
  public:
    std::vector<std::shared_ptr<Node<T>>>  list;
    Job *job;

    ForwardZddNode() : job(nullptr) {}

    ~ForwardZddNode() {
        list.clear();
    }


    std::shared_ptr<Node<T>> add_weight(int _weight, int _layer, bool _rootnode = false, bool _terminal_node = false){
        for (auto &it : list) {
            if (it->GetWeight() == _weight) {
                return it;
            }
        }

        std::shared_ptr<Node<T>> node = std::make_shared<Node<T>>(_weight, _layer, _rootnode, _terminal_node);
        node->set_job(job);
        list.push_back(std::move(node));
        return list.back();
    }

    void set_job(Job *_job){
        job = _job;
    }

    Job* get_job(){
        return job;
    }
};

template<typename T>
bool my_compare(const std::shared_ptr<Node<T>> &lhs, const std::shared_ptr<Node<T>> &rhs){
    return *lhs < *rhs;
}

template<typename E, typename T> class ForwardZddBase : public 
    tdzdd::DdEval<E, ForwardZddNode<T>, Optimal_Solution<T>> {
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

    void initializepi(T *_pi){
        pi = _pi;
    }

    virtual void initializenode(ForwardZddNode<T>& n) const  = 0;

    virtual void initializerootnode(ForwardZddNode<T>& n) const  = 0;

    virtual void evalNode(ForwardZddNode<T>& n) const = 0;

    // virtual Optimal_Solution<double> get_objective(ForwardZddNode<T>& n) const = 0;

    Optimal_Solution<T> get_objective(ForwardZddNode<T> &n) const {
        Optimal_Solution<T> sol(-pi[num_jobs]);

        int weight;


        auto m =  std::max_element(n.list.begin(), n.list.end(),my_compare<T>);
        weight = (*m)->get_weight();

        Label<T> *ptr_node = &((*m)->forward_label1);

        while(ptr_node->get_previous() != nullptr) {
            Label<T> *aux_prev_node = ptr_node->get_previous();
            Job *aux_job = aux_prev_node->get_job();
            if(ptr_node->get_high()) {
                sol.C_max += aux_job->processing_time;
                if(ptr_node->get_node()) {
                    sol.cost = sol.cost + value_Fj(ptr_node->get_weight(), aux_job) ;
                    sol.obj += pi[aux_job->job] - value_Fj(ptr_node->get_weight(), aux_job);
                } else {
                    sol.cost = sol.cost + value_Fj(weight, aux_job);
                    sol.obj += pi[aux_job->job] - value_Fj(weight, aux_job);
                }
                g_ptr_array_add(sol.jobs, aux_job);
            }
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

    void initializenode(ForwardZddNode<T>& n) const override {
        for (auto &it: n.list) {
            if(it->get_weight() == 0) {
                it->forward_label1.update_solution(-pi[num_jobs], nullptr, false);
                it->forward_label2.update_solution(-DBL_MAX/2, nullptr, false);
            } else {
                it->forward_label1.update_solution(-DBL_MAX/2, nullptr, false);
                it->forward_label2.update_solution(-DBL_MAX/2, nullptr, false);
            }
        }
    }

    void initializerootnode(ForwardZddNode<T> &n) const override {
        for(auto &it: n.list){
            it->forward_label1.f = -pi[num_jobs];
            it->forward_label2.set_f(-DBL_MAX/2);
        }
    }

    void evalNode(ForwardZddNode<T> &n) const override
    {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto &it : n.list) {
            int      weight = it->get_weight();
            T g;
            std::shared_ptr<Node<T>> p0 = it->n;
            std::shared_ptr<Node<T>> p1 = it->y;
            double result = - value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job];

            /**
             * High edge calculation
             */
            Job *prev = it->forward_label1.get_previous_job();
            Job *aux1 = p1->forward_label1.get_previous_job();
            double diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );

            if(prev != tmp_j && diff) {
                g = it->forward_label1.get_f() + result;
                if(g > p1->forward_label1.get_f()) {
                    if(aux1 != tmp_j) {
                        p1->forward_label2.update_solution(p1->forward_label1);
                    }
                    p1->forward_label1.update_solution(g, &(it->forward_label1), true);
                } else if ((g > p1->forward_label2.get_f()) && (aux1 != tmp_j)) {
                    p1->forward_label2.update_solution(g, &(it->forward_label1), true);
                }
            } else  {
                g = it->forward_label2.get_f() + result;
                prev = it->forward_label2.get_previous_job();
                diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );
                if(diff) {
                    if(g > p1->forward_label1.get_f()) {
                        if(aux1 != tmp_j) {
                            p1->forward_label2.update_solution(p1->forward_label1);
                        }
                        p1->forward_label1.update_solution(g, &(it->forward_label2), true);
                    } else if ((g > p1->forward_label2.get_f()) && (aux1 != tmp_j)) {
                        p1->forward_label2.update_solution(g, &(it->forward_label2), true);
                    }
                }
            }

            /**
             * Low edge calculation
             */
            aux1 = p0->forward_label1.get_previous_job();
            if(it->forward_label1.get_f() > p0->forward_label1.get_f()) {
                if(prev != aux1) {
                    p0->forward_label2.update_solution(p0->forward_label1);
                }
                p0->forward_label1.update_solution(it->forward_label1);
                if(it->forward_label2.get_f() > p0->forward_label2.get_f()) {
                    p0->forward_label2.update_solution(it->forward_label2);
                }
            } else if ((it->forward_label1.get_f() > p0->forward_label2.get_f()) && (aux1 != prev)){
                p0->forward_label2.update_solution(it->forward_label1);
            } else if ((it->forward_label2.get_f() > p0->forward_label2.get_f())) {
                p0->forward_label2.update_solution(it->forward_label2);
            }
        }
    }

    Optimal_Solution<T> getValue(ForwardZddNode<T> const &n) {
        Optimal_Solution<T> sol;

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

    void initializenode(ForwardZddNode<T>& n) const override {
        for (auto &it: n.list) {
            if(it->get_weight() == 0) {
                it->forward_label1.update_solution(-pi[num_jobs], nullptr, false);
            } else {
                it->forward_label1.update_solution(-DBL_MAX/2, nullptr, false);
            }
        }
    }

    void initializerootnode(ForwardZddNode<T> &n) const override {
        for(auto &it: n.list){
            it->forward_label1.f = -pi[num_jobs];
        }
    }

    void initializepi(T *_pi){
        pi = _pi;
    }

    void evalNode(ForwardZddNode<T> &n) const override {
        Job *tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto &it : n.list) {
            int      weight = it->get_weight();
            T g;
            std::shared_ptr<Node<T>> p0 = it->n;
            std::shared_ptr<Node<T>> p1 = it->y;
            double result = - value_Fj(weight + tmp_j->processing_time, tmp_j) + pi[tmp_j->job];

            /**
             * High edge calculation
             */
            g = it->forward_label1.get_f() + result;
            if(g > p1->forward_label1.get_f()) {
                p1->forward_label1.update_solution(g, &(it->forward_label1), true);
            }

            /**
             * Low edge calculation
             */
            if(it->forward_label1.get_f() > p0->forward_label1.get_f()) {
                p0->forward_label1.update_solution(it->forward_label1);
            }
        }
    }

    Optimal_Solution<T> getValue(ForwardZddNode<T> const &n) {
        Optimal_Solution<T> sol;
        return sol;
    }
};

#endif // DURATION_ZDD_HPP