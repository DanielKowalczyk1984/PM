#ifndef DURATION_ZDD_HPP
#define DURATION_ZDD_HPP
#include <tdzdd/DdEval.hpp>

#include <OptimalSolution.hpp>
#include <node_duration.hpp>

#include <vector>
#include <algorithm>
using namespace std;

template<typename T>
class PricerDurationZDD {
  public:

    PricerDurationZDD() {}

    ~PricerDurationZDD() {
        list.clear();
    }

    std::vector<NodeDuration<T>> list;

    NodeDuration<T>* add_weight(int _weight, int _layer, bool _rootnode = false, bool _terminal_node = false){
        for (NodeDuration<T> &it : list) {
            if (it.GetWeight() == _weight) {
                return &it;
            }
        }

        NodeDuration<T> node(_weight, _layer, _rootnode, _terminal_node);
        list.push_back(node);
        return &(list[list.size()]);
    }
};

// template<typename T>
// bool compare_func(const NodeDuration<T> &l, const NodeDuration<T> &r){
//     return l < r;
// }

template<typename E, typename T> class DurationZDD : public
    tdzdd::DdEval<E, PricerDurationZDD<T>, Optimal_Solution<T>> {
  private:
    T *pi;
    GPtrArray *interval_list;
    int num_layers;
    int num_jobs;
  public:
    DurationZDD(T *_pi, GPtrArray *_interval_list, int _num_jobs): pi(_pi),
        interval_list(_interval_list), num_jobs(_num_jobs) {
        num_layers = interval_list->len;
    }

    void initializenode(PricerDurationZDD<T>& n){
        for (auto &it: n.list) {
            it.prev1.SetF(-DBL_MAX);
            it.prev2.SetF(-DBL_MAX);
            it.prev1.SetPrev(nullptr);
            it.prev2.SetPrev(nullptr);
            it.prev1.SetHigh(false);
            it.prev2.SetHigh(false);
        }
    }

    void initializerootnode(PricerDurationZDD<T> &n) {
        for(auto &it: n.list){
            it.prev1.SetF(pi[num_jobs]);
            it.prev2.SetF(-DBL_MAX);
        }
    }

    void evalNode(PricerDurationZDD<T> &n, int i) const
    {
        int j = num_layers - i;
        assert(j >= 0 && j <= num_layers - 1);
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          interval_list, j);
        Job *tmp_j = tmp_pair->j;
        double result;

        for (NodeDuration<T>& it : n.list) {
            int      weight = it.GetWeight();
            T g;
            NodeDuration<T>* p0 = it.n;
            NodeDuration<T>* p1 = it.y;
            result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

            /**
             * High edge calculation
             */
            if(it.prev1.GetJob() != tmp_j) {
                g = it.prev1.GetF() + result;
                if(g > p1->prev1.GetF()) {
                    if(p1->prev1.GetJob() != tmp_j) {
                        p1->prev2.SetF(p1->prev1.GetF());
                        p1->prev2.SetPrev(p1->prev1.GetPrev());
                        p1->prev2.SetHigh(p1->prev1.GetHigh());
                    }
                    p1->prev1.SetF(g);
                    p1->prev1.SetPrev(&(it.prev1));
                    p1->prev1.SetHigh(true);
                } else if ((g > p1->prev2.GetF()) && (p1->prev2.GetJob() != tmp_j)) {
                    p1->prev2.SetF(g);
                    p1->prev2.SetPrev(&(it.prev1));
                    p1->prev2.SetHigh(true);
                }

            } else {
                g = it.prev2.GetF() + result;
                if(g > p1->prev1.GetF()) {
                    if(p1->prev1.GetJob() != tmp_j) {
                        p1->prev2.SetF(p1->prev1.GetF());
                        p1->prev2.SetPrev(p1->prev1.GetPrev());
                        p1->prev2.SetHigh(p1->prev1.GetHigh());
                    }
                    p1->prev1.SetF(g);
                    p1->prev1.SetPrev(&(it.prev2));
                    p1->prev1.SetHigh(true);
                } else if ((g > p1->prev2.GetF()) && (p1->prev2.GetJob() != tmp_j)) {
                    p1->prev2.SetF(g);
                    p1->prev2.SetPrev(&(it.prev2));
                    p1->prev2.SetHigh(true);
                }
            }

            /**
             * Low edge calculation
             */
            if(it.prev1.GetF() > p0->prev1.GetF()) {
                if(it.prev1.GetJob() != p0->prev1.GetJob()) {
                    p0->prev2.SetF(p0->prev1.GetF());
                    p0->prev2.SetPrev(p0->prev1.GetPrev());
                    p0->prev2.SetHigh(p0->prev1.GetHigh());
                }
                p0->prev1.SetF(it.prev1.GetF());
                p0->prev1.SetPrev(it.prev1.GetPrev());
                p0->prev1.SetHigh(it.prev1.GetHigh());
            } else if ((it.prev1.GetF() > p0->prev2.GetF()) && (it.prev1.GetJob() != p0->prev1.GetJob())){
                p0->prev2.SetF(it.prev1.GetF());
                p0->prev2.SetPrev(it.prev1.GetPrev());
                p0->prev2.SetHigh(it.prev1.GetHigh());
            }
        }
    }


    Optimal_Solution<T> get_objective(PricerDurationZDD<T> *n) {
        Optimal_Solution<T> sol;
        Job *aux;
        sol.cost = 0;
        sol.C_max = 0;
        typename std::vector<NodeDuration<T>>::iterator m =  std::max_element(n->list.begin(), n->list.end());
        max
        return sol;
    }
};

#endif // DURATION_ZDD_HPP