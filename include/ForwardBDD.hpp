#ifndef FORWARD_BDD_HPP
#define FORWARD_BDD_HPP
#include <tdzdd/DdEval.hpp>
#include <OptimalSolution.hpp>
#include <node_duration.hpp>


template<typename E, typename T> class ForwardBddBase : public 
    tdzdd::DdEval<E, Node<T>, Optimal_Solution<T>> {
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

    virtual void initializenode(Node<T>& n) const  = 0;

    virtual void initializerootnode(Node<T>& n) const  = 0;

    virtual void evalNode(Node<T>& n) const = 0;

    Optimal_Solution<T> get_objective(Node<T> &n) const {
        Optimal_Solution<T> sol(-pi[num_jobs]);

        PrevNode<T> *ptr_node = &(n.prev1);

        while(ptr_node->GetPrev() != nullptr) {
            PrevNode<T> *aux_prev_node = ptr_node->GetPrev();
            Job *aux_job = aux_prev_node->GetJob();
            if(ptr_node->GetHigh()) {
                sol.C_max += aux_job->processingime;
                sol.push_job_back(aux_job, aux_prev_node->GetWeight(), pi[aux_job->job]);
            }
            ptr_node = aux_prev_node;
        }

        // weight = 0;
        // for(int i = sol.jobs->len - 1; i >= 0; i--){
        //     aux_job = (Job *) g_ptr_array_index(sol.jobs, i);
        //     weight += aux_job->processingime;
        //     sol.cost +=  value_Fj(weight, aux_job);
        //     sol.obj += pi[aux_job->job] - value_Fj(weight, aux_job);
        // }
        // printf("test weight 2 = %d\n", weight);

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
        if(n.GetWeight() == 0) {
            n.prev1.UpdateSolution(-pi[num_jobs], nullptr, false);
            n.prev2.UpdateSolution(-DBL_MAX/2, nullptr, false);
        } else {
            n.prev1.UpdateSolution(-DBL_MAX/2, nullptr, false);
            n.prev2.UpdateSolution(-DBL_MAX/2, nullptr, false);
        }
    }

    void initializerootnode(Node<T> &n) const override {
        n.prev1.f = pi[num_jobs];
        n.prev2.SetF(-DBL_MAX/2);
    }

    void evalNode(Node<T> &n) const override
    {
        Job *tmp_j = n.GetJob();
        assert(tmp_j != nullptr);
        double result;
        bool diff;

        int      weight = n.GetWeight();
        T g;
        Node<T>* p0 = n.child[0];
        Node<T>* p1 = n.child[1];
        result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

        /**
         * High edge calculation
         */
        Job *prev = n.prev1.GetPrevJob();
        Job *aux1 = p1->prev1.GetPrevJob();
        diff = (prev == nullptr ) ? true : (value_diff_Fij(weight, tmp_j, prev) >= 0 );

        if(prev != tmp_j && diff) {
            g = n.prev1.GetF() + result;
            if(g > p1->prev1.GetF()) {
                if(aux1 != tmp_j) {
                    p1->prev2.UpdateSolution(p1->prev1);
                }
                p1->prev1.UpdateSolution(g, &(n.prev1), true);
            } else if ((g > p1->prev2.GetF()) && (aux1 != tmp_j)) {
                p1->prev2.UpdateSolution(g, &(n.prev1), true);
            }
        } else  {
            g = n.prev2.GetF() + result;
            if(g > p1->prev1.GetF()) {
                if(aux1 != tmp_j) {
                    p1->prev2.UpdateSolution(p1->prev1);
                }
                p1->prev1.UpdateSolution(g, &(n.prev2), true);
            } else if ((g >= p1->prev2.GetF()) && (aux1 != tmp_j)) {
                p1->prev2.UpdateSolution(g, &(n.prev2), true);
            }
        }

        /**
         * Low edge calculation
         */
        aux1 = p0->prev1.GetPrevJob();
        if(n.prev1.GetF() > p0->prev1.GetF()) {
            if(prev != aux1) {
                p0->prev2.UpdateSolution(p0->prev1);
            }
            p0->prev1.UpdateSolution(n.prev1);
        } else if ((n.prev1.GetF() > p0->prev2.GetF()) && (aux1 != prev)){
            p0->prev2.UpdateSolution(n.prev1);
        }
    }

    Optimal_Solution<T> getValue(Node<T> const &n){
        Optimal_Solution<T> sol;

        return sol;
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
        if(n.GetWeight() == 0) {
            n.prev1.UpdateSolution(-pi[num_jobs], nullptr, false);
        } else {
            n.prev1.UpdateSolution(-DBL_MAX/2, nullptr, false);
        }
    }

    void initializerootnode(Node<T> &n) const override {
        n.prev1.f = pi[num_jobs];
    }

    void initializepi(T *_pi){
        pi = _pi;
    }

    void evalNode(Node<T> &n) const override {
        Job *tmp_j = n.GetJob();
        assert(tmp_j != nullptr);
        T result;

        int      weight = n.GetWeight();
        T g;
        Node<T>* p0 = n.child[0];
        Node<T>* p1 = n.child[1];
        result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

        /**
         * High edge calculation
         */
        g = n.prev1.GetF() + result;
        if(g > p1->prev1.GetF()) {
            p1->prev1.UpdateSolution(g, &(n.prev1), true);
        }

        /**
         * Low edge calculation
         */
        if(n.prev1.GetF() > p0->prev1.GetF()) {
            p0->prev1.UpdateSolution(n.prev1);
        }
    }

    Optimal_Solution<T> getValue(Node<T> const &n){
        Optimal_Solution<T> sol;
        return sol;
    }
};

#endif // FORWARD_BDD_HPP
