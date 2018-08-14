#ifndef WEIGHT_ZDD_HPP
#define WEIGHT_ZDD_HPP
#include <tdzdd/DdEval.hpp>
#include <tdzdd/dd/NodeTable.hpp>

#include <node.hpp>
#include <OptimalSolution.hpp>

template <typename T>
class PricerWeightZDD {
public:
    std::vector<std::shared_ptr<node<T>>> list;

    PricerWeightZDD() {};

    ~PricerWeightZDD() noexcept{
        list.clear();
    }

    std::shared_ptr<node<T>>& add_weight(int _weight, int layer, int njobs)
    {
        for (auto& it : list) {
            if (it->weight == _weight || it->terminal) {
                return it;
            }
        }

        std::shared_ptr<node<T>> p = std::make_shared<node<T>>();
        p->layer = layer;
        p->weight = _weight;
        list.push_back(p);
        return (list.back());
    }

    void add_terminal_node(int j, int layer, int njobs){
        std::shared_ptr<node<T>> p = std::make_shared<node<T>>();
        if(j == 0){
            p->obj = -DBL_MAX;
        } else {
            p->obj = -DBL_MAX;
        }

        p->take = false;
        p->terminal = true;
        p->layer = layer;
        list.push_back(p);
    }

    void init_terminal_node(int j, int njobs)
    {
        for (auto& it : list) {
            it->obj = -DBL_MAX;
            it->take = false;
            it->prev = nullptr;
        }
    }
};


template <typename E, typename T>
class WeightZDD
    : public tdzdd::DdEval<E, PricerWeightZDD<T>, Optimal_Solution<T>> {
    T   *pi;
    GPtrArray *interval_list;
    int nlayers;
    int nbjobs;

public:
    WeightZDD(T *_pi, GPtrArray *_interval_list, int _nbjobs): pi(_pi),
        interval_list(_interval_list), nbjobs(_nbjobs)
    {
        nlayers = (int) interval_list->len;
    };

    void evalTerminal(PricerWeightZDD<T>& n)
    {
        for (auto &it : n.list) {
            it->obj = pi[nbjobs];
            it->c = pi[nbjobs];
            it->b = pi[nbjobs];
            it->prev_job = nullptr;
            it->prev_node = nullptr;
        }
    }

    void evalNode(PricerWeightZDD<T> *n, int i,
                  tdzdd::DdValues<PricerWeightZDD<T>, 2>& values) const
    {
        int j = nlayers - i;
        assert(j >= 0 && j <= nlayers - 1);
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          interval_list, j);
        Job *tmp_j = tmp_pair->j;
        double result;

        for (auto& it : n->list) {
            int      weight = it->weight;
            std::shared_ptr<node<T>> p0 = it->n;
            std::shared_ptr<node<T>> p1 = it->y;
            std::shared_ptr<node<T>> aux;
            T        obj0 = p0->obj;
            T        obj1 = p1->obj;
            result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

            it->obj1 = obj1 + result;
            it->obj0 = obj0;

            if ((it->calc)) {
                if (it->obj0 < it->obj1) {
                    it->obj = it->obj1;
                    it->take = true;
                    it->prev_job = tmp_j;
                } else {
                    it->obj = it->obj0;
                    it->take = false;
                    it->prev_job = p0->prev_job;
                }

                if (it->b < it->obj1) {
                    it->b = it->obj1;
                }

                if(it->c < it->obj0) {
                    it->c = it->obj0;
                }
            } else {
                it->obj = it->obj0;
                it->take = false;
                it->prev_job = p0->prev_job;
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
        }
    }

    Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram,
                                      tdzdd::DataTable<PricerWeightZDD<T>> *data_table,
                                      const tdzdd::NodeId *f)
    {
        Optimal_Solution<T> sol;
        sol.obj = (*data_table)[f->row()][f->col()].list[sol.C_max]->obj;
        std::shared_ptr<node<T>> ptr = (*data_table)[f->row()][f->col()].list[sol.C_max];
        job_interval_pair *tmp_pair;
        Job *tmp_j;


        while (ptr->layer != nlayers) {
            tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, ptr->layer);
            tmp_j = tmp_pair->j;

            if (ptr->obj == ptr->obj1) {
                g_ptr_array_add(sol.jobs, tmp_j);
                // auto e = ptr->out_edge[0].lock();
                // g_ptr_array_add(sol.e_list, &(e->id));
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                ptr = ptr->y;
            } else {
                // auto e = ptr->out_edge[1].lock();
                // g_ptr_array_add(sol.e_list, &(e->id));
                ptr = ptr->n;
            }
        }

        return sol;
    }
};


#endif // WEIGHT_ZDD_HPP


