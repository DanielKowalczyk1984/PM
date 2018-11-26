#ifndef DURATION_BDD_HPP
#define DURATION_BDD_HPP
#include <tdzdd/DdEval.hpp>
#include <cfloat>
#include <vector>
#include <iostream>
#include <StateNode.hpp>
#include <OptimalSolution.hpp>


template<typename T>
class PricerInfoBDD {
  public:
    T obj;
    std::vector<StateNode<T>> bucket;
    std::vector<Job*> jobs;
    int sum_w;
    int sum_p;
    int cost;

    PricerInfoBDD &operator=(const PricerInfoBDD &other) {
        sum_w = other.sum_w;
        sum_p = other.sum_p;
        cost = other.cost;
        obj = other.obj;
        return *this;
    };

    friend std::ostream &operator<<(std::ostream &os, PricerInfoBDD<T> const &o) {
        os << "max = " << o.obj << "," << std::endl << "cost = " << o.cost << std::endl;
        return os;
    };

    void init_terminal_node(int one) {
        obj = one ? -DBL_MAX : -DBL_MAX;
        jobs.reserve(40);
        sum_p = 0;
    }

    void init_node(int weight, int njobs) {
        obj = -DBL_MAX;
        jobs.reserve(njobs);
        sum_p = weight;
    };
};


template<typename E, typename T>
class DurationBDD: public
    tdzdd::DdEval<E, PricerInfoBDD<T>, Optimal_Solution<T> > {
    T *pi;
    GPtrArray *interval_list;
    int nlayers;
    int nbjobs;


  public:
    DurationBDD(T *_pi, GPtrArray *_interval_list, int _nbjobs)
        : pi(_pi), interval_list(_interval_list), nbjobs(_nbjobs) {
            nlayers = interval_list->len;
    };

    // void evalTerminal(PricerInfoBDD<T> &n, bool one) {
    //     n.obj = one ? 0 : -1871286761.0;
    //     n.cost = 0;
    //     n.sum_p = 0;
    //     n.x.resize(0);
    // }

    void evalNode(PricerInfoBDD<T> &n) const {
        int j = nlayers ;
        assert(j >= 0 && j <= nlayers - 1);
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          interval_list, j);
        Job *tmp_j = tmp_pair->j;
        PricerInfoBDD<T> *n0 ;
        PricerInfoBDD<T> *n1 ;


        if (n0->obj < n.obj) {
            n0->obj = n.obj;
            n0->cost = n.cost;
            n0->sum_p = n.sum_p;
            n0->jobs = n.jobs;
        }

        if (n1->obj < n.obj - (T) value_Fj(n.sum_p + tmp_j->processingime, tmp_j) +  pi[tmp_j->job] ) {
            n1->obj = n.obj - (T) value_Fj(n.sum_p + tmp_j->processingime, tmp_j) +  pi[tmp_j->job];
            n1->cost = n.cost + value_Fj(n.sum_p + tmp_j->processingime, tmp_j);
            n1->sum_p = n.sum_p + tmp_j->processingime;
            n1->jobs = n.jobs;
            n1->jobs.push_back(tmp_j);
        }
    }

    void initializenode(PricerInfoBDD<T> &n) {
        n.obj = -DBL_MAX;
        n.cost = 0;
        n.jobs.clear();
        n.jobs.reserve(nbjobs);
    }

    void initializerootnode(PricerInfoBDD<T> &n) {
        n.obj = pi[nbjobs];
        n.cost = 0;
        n.jobs.clear();
        n.jobs.reserve(nbjobs);
    }

    Optimal_Solution<T> get_objective(PricerInfoBDD<T> &n) {
        Optimal_Solution<T> sol;
        sol.cost = 0;
        sol.obj = n.obj;
        sol.C_max = 0;
        sol.jobs = g_ptr_array_sized_new(n.jobs.size());
        for(const auto& it: n.jobs){
            g_ptr_array_add(sol.jobs, it);
            sol.C_max += ((Job *) it)->processingime;
            sol.cost += value_Fj(sol.C_max, it);
        }
        return sol;
    }


};


#endif // DURATION_BDD_HPP


