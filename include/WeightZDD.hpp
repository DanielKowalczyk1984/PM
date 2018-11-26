#ifndef WEIGHT_ZDD_HPP
#define WEIGHT_ZDD_HPP
#include <tdzdd/DdEval.hpp>
#include <tdzdd/dd/NodeTable.hpp>

#include <node.hpp>
#include <OptimalSolution.hpp>

template <typename T> class PricerWeightZDD {
public:
    std::vector<std::shared_ptr<node<T>>> list;
    Job *job;

    PricerWeightZDD() {};

    ~PricerWeightZDD() noexcept{
        list.clear();
    }

    std::shared_ptr<node<T>>& add_weight(int _weight, int layer)
    {
        for (auto& it : list) {
            if (it->weight == _weight || it->terminal) {
                return it;
            }
        }

        std::shared_ptr<node<T>> p{std::make_shared<node<T>>()};
        p->job = job;
        p->weight = _weight;
        list.push_back(std::move(p));
        return (list.back());
    }

    void set_job(Job *_job){
        job = _job;
    }

    Job* get_job(){
        return job;
    }
};


template <typename E, typename T>
class WeightZDD : public tdzdd::DdEval<E, PricerWeightZDD<T>, Optimal_Solution<T>> {
    int nbjobs;
    T   *pi;

public:

    WeightZDD() : nbjobs(0), pi(nullptr){
    }

    WeightZDD(T *_pi,  int _nbjobs):  nbjobs(_nbjobs), pi(_pi){
    };

    WeightZDD(int _nbjobs): nbjobs(_nbjobs){
    };

    WeightZDD(const WeightZDD<E, T> &src){
        nbjobs = src.nbjobs;
        pi = src.pi;
    }

    void evalTerminal(PricerWeightZDD<T>& n)
    {
        for (auto &it : n.list) {
            it->obj = pi[nbjobs];
            it->c = pi[nbjobs];
            it->b = pi[nbjobs];
            it->terminal = true;
        }
    }

    void initializepi(T *_pi){
        pi = _pi;
    }

    void evalNode(PricerWeightZDD<T> *n) const
    {
        Job *tmp_j = n->get_job();
        assert(tmp_j != nullptr);

        T result, aux_obj0, aux_obj1;

        for (auto& it : n->list) {
            int      weight = it->weight;
            std::shared_ptr<node<T>> p0 = it->n;
            std::shared_ptr<node<T>> p1 = it->y;
            T        obj0 = p0->obj;
            T        obj1 = p1->obj;
            result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

            aux_obj1 = obj1 + result;
            aux_obj0 = obj0;

            if ((it->calc)) {
                if (aux_obj0 < aux_obj1) {
                    it->obj = aux_obj1;
                    it->take = true;
                } else {
                    it->obj = aux_obj0;
                    it->take = false;
                }

                if (it->b < aux_obj1) {
                    it->b = aux_obj1;
                }

                if(it->c < aux_obj0) {
                    it->c = aux_obj0;
                }
            } else {
                it->obj = aux_obj0;
                it->take = false;
            }
        }
    }

    void initializenode(PricerWeightZDD<T>& n)
    {
        for (auto& it : n.list) {
            it->obj = -DBL_MAX;
            it->b = pi[nbjobs];
            it->c = pi[nbjobs];
            it->dist = -DBL_MAX;
            it->take = false;
            it->terminal = false;
        }
    }

    Optimal_Solution<T> get_objective(PricerWeightZDD<T> &n)
    {
        Optimal_Solution<T> sol;
        std::shared_ptr<node<T>> ptr = n.list.front();
        sol.obj = pi[nbjobs];

        do {
            Job *tmp_j = ptr->job;
            if (ptr->take) {
                g_ptr_array_add(sol.jobs, tmp_j);
                // auto e = ptr->out_edge[0].lock();
                // g_ptr_array_add(sol.e_list, &(e->id));
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                sol.obj += -value_Fj(sol.C_max, tmp_j) + pi[tmp_j->job];
                ptr = ptr->y;
            } else {
                // auto e = ptr->out_edge[1].lock();
                // g_ptr_array_add(sol.e_list, &(e->id));
                ptr = ptr->n;
            }
        }while (ptr->job);

        return sol;
    }
};


#endif // WEIGHT_ZDD_HPP


