#include <solution.h>
#include <algorithm>
#include <array>
#include <boost/dynamic_bitset.hpp>
#include <cfloat>
#include <climits>
#include <unordered_map>
#include <vector>
#include "tdzdd/DdEval.hpp"
#include "tdzdd/dd/NodeTable.hpp"
#include <solution.h>
#include <interval.h>

template <typename T>
class node {
   public:
    int      layer;
    int      weight;
    T        obj;
    bool     take;
    node<T> *y;
    node<T> *n;

    node()
        : layer(0),
          weight(0),
          obj(0.0),
          take(false),
          y(nullptr),
          n(nullptr){

          };

    friend std::ostream &operator<<(std::ostream &os, node<T> const &o) {
        os << "layer = " << o.layer << ", weight = " << o.weight
           << ", obj = " << o.obj << std::endl;

        return os;
    };

    node<T>(const node<T> &other) {
        layer = other.layer;
        weight = other.weight;
        obj = other.obj;
        take = other.take;
        y = other.y;
        n = other.n;
    };
};

// /**
//  * BDD
//  */
template <typename T>
class PricerWeightBDD {
   public:
    T    obj;
    bool take;
    int  sum_w;
    int  sum_p;
    int  cost;

    PricerWeightBDD() : obj(0), take(false), sum_w(0), sum_p(0), cost(0) {}

    PricerWeightBDD &operator=(const PricerWeightBDD &other) {
        sum_w = other.sum_w;
        sum_p = other.sum_p;
        cost = other.cost;
        obj = other.obj;
        take = other.take;
        return *this;
    };

    friend std::ostream &operator<<(std::ostream &            os,
                                    PricerWeightBDD<T> const &o) {
        os << "max = " << o.obj << "," << std::endl
           << "cost = " << o.cost << std::endl;
        return os;
    }

    void init_terminal_node(int one) {
        obj = one ? 0.0 : -1871286761.0;
        sum_p = 0;
        take = false;
    }

    void init_node(int weight) {
        obj = 0.0;
        sum_p = weight;
        take = false;
    };
};

template <typename T>
using my_iterator = typename std::vector<node<T> *>::iterator;

template <typename T>
class PricerWeightZDD {
   public:
    std::vector<node<T> *> list;

    PricerWeightZDD(){};

    ~PricerWeightZDD() {
        for (my_iterator<T> it = list.begin(); it != list.end(); it++) {
            delete *it;
        }
    }

    node<T> *add_weight(int _weight, int layer) {
        for (my_iterator<T> it = list.begin(); it != list.end(); it++) {
            if ((*it)->weight == _weight) {
                return (*it);
            }
        }

        node<T> *p = new node<T>();
        p->layer = layer;
        p->weight = _weight;
        p->obj = 0.0;
        p->take = false;
        list.push_back(p);

        return (list.back());
    }

    void init_terminal_node(int j) {
        for (my_iterator<T> it = list.begin(); it != list.end(); it++) {
            (*it)->obj = j ? 0.0 : -1871286761.0;
            (*it)->take = false;
        }
    }
};

template <typename T>
class PricerFarkasZDD {
   public:
    T    obj;
    bool take;

    PricerFarkasZDD() : obj(0), take(0){};

    ~PricerFarkasZDD(){};

    void init_terminal_node(int one) { obj = one ? 0.0 : 1871286761.0; }

    void init_node() { take = false; }
};

template <typename T>
class Optimal_Solution {
   public:
    T                obj;
    int              cost;
    int              C_max;
    std::vector<int> jobs;
    // Optimal_Solution(PricerInfoZDD<T> *node, const int L) {
    //     int                      it_max = node->get_max(L);
    //     boost::dynamic_bitset<> *ptr_max = &node->A[it_max];
    //     size_t                   it = ptr_max->find_first();
    //     obj = node->obj[it_max];
    //     cost = node->cost[it_max];
    //     C_max = 0;

    //     while (it != boost::dynamic_bitset<>::npos) {
    //         if ((*ptr_max)[it]) {
    //             jobs.push_back(it);
    //         }

    //         it = ptr_max->find_next(it);
    //     }
    // }

    // Optimal_Solution(PricerInfoBDD<T> *node) : jobs(node->jobs) {
    //     obj = node->obj;
    //     // jobs = std::vector<int>(node->jobs.begin(),node->jobs.end());
    //     cost = node->cost;
    //     C_max = 0;
    // }

    Optimal_Solution() : obj(0), cost(0), C_max(0) {}

    Optimal_Solution &operator=(const Optimal_Solution &other) {
        obj = other.obj;
        cost = other.cost;
        jobs = other.jobs;
        C_max = other.C_max;
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &             os,
                                    Optimal_Solution<T> const &o) {
        os << "obj = " << o.obj << "," << std::endl
           << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;
        std::vector<int>::const_iterator it = o.jobs.begin();

        for (; it != o.jobs.end(); ++it) {
            std::cout << *it << " ";
        }

        std::cout << std::endl;
        return os;
    };
};

template <typename E, typename T>
class WeightBDD
    : public tdzdd::DdEval<E, PricerWeightBDD<T>, Optimal_Solution<T>> {
    T *  pi;
    GPtrArray *interval_list;
    int nlayers;
    int nbjobs;
    job_interval_pair *tmp_pair;
    Job *tmp_j;
    interval *tmp_interval;

   public:
    WeightBDD(T *_pi, GPtrArray *_interval_list, int _nbjobs)
        : pi(_pi), interval_list(_interval_list), nbjobs(_nbjobs){};

    void evalTerminal(PricerWeightBDD<T> &n) {
        n.obj = pi[nbjobs];
        n.cost = 0.0;
        n.sum_w = 0.0;
        n.sum_p = 0.0;
        n.take = false;
    }

    void evalNode(PricerWeightBDD<T> *n,
                  int                 i,
                  tdzdd::DdValues<PricerWeightBDD<T>, 2> &values) const {
        int j = nlayers - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerWeightBDD<T> *n0 = values.get_ptr(0);
        PricerWeightBDD<T> *n1 = values.get_ptr(1);
        tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
        tmp_j = tmp_pair->j;

        if (n0->obj >= n1->obj - value_Fj(n->sum_p + tmp_j->processingime, tmp_j) + pi[tmp_j->job]) {
            n->obj = n0->obj;
            n->take = false;
        } else {
            n->obj = n1->obj - value_Fj(n->sum_p + tmp_j->processingime, tmp_j) + pi[tmp_j->job];
            n->take = true;
        }
    }

    void initializenode(PricerWeightBDD<T> &n) { n.take = false; }

    Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram,tdzdd::DataTable<PricerWeightBDD<T>> *data_table,
        const tdzdd::NodeId *f) {
        Optimal_Solution<T> sol;
        sol.obj = (*data_table)[f->row()][f->col()].obj;
        sol.cost = 0;
        sol.C_max = 0;
        tdzdd::NodeId cur_node = *f;
        int j = nlayers - cur_node.row();
        tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
        tmp_j = tmp_pair->j;

        while (cur_node.row() != 0 || cur_node.col() != 0) {
            if ((*data_table)[cur_node.row()][cur_node.col()].take) {
                sol.jobs.push_back(tmp_j->job);
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                cur_node = diagram.privateEntity().child(cur_node, 1);
            } else {
                cur_node = diagram.privateEntity().child(cur_node, 0);
            }

            j = nbjobs - cur_node.row();
            tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
            tmp_j = tmp_pair->j;
        }

        return sol;
    }
};

template <typename E, typename T>
class WeightZDD
    : public tdzdd::DdEval<E, PricerWeightZDD<T>, Optimal_Solution<T>> {
    T *  pi;
    GPtrArray *interval_list;
    int nlayers;
    int nbjobs;
    job_interval_pair *tmp_pair;
    Job *tmp_j;
    interval *tmp_interval;

   public:
    WeightZDD(T *_pi, GPtrArray *_interval_list, int _nbjobs): pi(_pi), interval_list(_interval_list), nbjobs(_nbjobs) {
            nlayers = (int) interval_list->len;
    };

    void evalTerminal(PricerWeightZDD<T> &n) {
        for (my_iterator<T> it = n.list.begin(); it != n.list.end(); it++) {
            (*it)->obj = pi[nbjobs];
        }
    }

    void evalNode(PricerWeightZDD<T> *n,int i, tdzdd::DdValues<PricerWeightZDD<T>, 2> &values) const {
        int j = nlayers - i;
        assert(j >= 0 && j <= nlayers - 1);
        tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
        tmp_j = tmp_pair->j;

        for (my_iterator<T> it = n->list.begin(); it != n->list.end(); it++) {
            int      weight = ((*it)->weight);
            node<T> *p0 = (*it)->n;
            node<T> *p1 = (*it)->y;
            T        obj0 = p0->obj;
            T        obj1 = p1->obj;

            if (obj0 >= obj1 - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job]) {
                (*it)->obj = obj0;
                (*it)->take = false;
            } else {
                (*it)->obj = obj1 - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];
                (*it)->take = true;
            };
        }
    }

    void initializenode(PricerWeightZDD<T> &n) {
        for (my_iterator<T> it = n.list.begin(); it != n.list.end(); it++) {
            (*it)->obj = 0.0;
            (*it)->take = false;
        }
    }

    Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram,tdzdd::DataTable<PricerWeightZDD<T>> *data_table,
        const tdzdd::NodeId *f) {
        Optimal_Solution<T> sol;
        sol.C_max = 0;
        sol.cost = 0;
        sol.obj = (*data_table)[f->row()][f->col()].list[sol.C_max]->obj;
        node<T> *ptr = (*data_table)[f->row()][f->col()].list[sol.C_max];
        int j = ptr->layer;
        job_interval_pair *tmp_pair;
        Job *tmp_j;

        while (ptr->layer != nlayers) {
            tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
            tmp_j = tmp_pair->j;
            if (ptr->take) {
                sol.jobs.push_back(tmp_j->job);
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                ptr = ptr->y;
            } else {
                ptr = ptr->n;
            }
            j = ptr->layer;
        }

        return sol;
    }
};

/**
 * Farkas
 */
template <typename E, typename T>
class FarkasZDD
    : public tdzdd::DdEval<E, PricerFarkasZDD<T>, Optimal_Solution<T>> {
    T *  pi;
    Job *jobarray;
    int  nbjobs;
    int  H_min;
    int  H_max;

   public:
    FarkasZDD(T *_pi, Job *_jobarray, int _nbjobs, int Hmin, int Hmax)
        : pi(_pi),
          jobarray(_jobarray),
          nbjobs(_nbjobs),
          H_min(Hmin),
          H_max(Hmax){};

    void evalTerminal(PricerFarkasZDD<T> &n) { n.obj = pi[nbjobs]; }

    void evalNode(PricerFarkasZDD<T> *n,
                  int                 i,
                  tdzdd::DdValues<PricerFarkasZDD<T>, 2> &values) const {
        int j = nbjobs - i;
        assert(j >= 0 && j <= nbjobs - 1);
        PricerFarkasZDD<T> *n0 = values.get_ptr(0);
        PricerFarkasZDD<T> *n1 = values.get_ptr(1);

        if (n0->obj < n1->obj + pi[j]) {
            n->obj = n0->obj;
            n->take = false;
        } else {
            n->obj = n1->obj + pi[j];
            n->take = true;
        }
    }

    void initializenode(PricerFarkasZDD<T> &n) { n.take = false; }

    Optimal_Solution<T> get_objective(
        tdzdd::NodeTableHandler<2>            diagram,
        tdzdd::DataTable<PricerFarkasZDD<T>> *data_table,
        const tdzdd::NodeId *                 f) {
        Optimal_Solution<T> sol;
        sol.obj = (*data_table)[f->row()][f->col()].obj;
        sol.cost = 0;
        sol.C_max = 0;
        tdzdd::NodeId cur_node = *f;
        int           j = nbjobs - cur_node.row();

        while (cur_node.row() != 0) {
            if ((*data_table)[cur_node.row()][cur_node.col()].take &&
                jobarray[j].releasetime <= sol.C_max &&
                sol.C_max + jobarray[j].processingime <= jobarray[j].duetime) {
                sol.jobs.push_back(j);
                sol.C_max += jobarray[j].processingime;
                sol.cost += jobarray[j].weight * sol.C_max;
                cur_node = diagram.privateEntity().child(cur_node, 1);
                j = nbjobs - cur_node.row();
            } else {
                cur_node = diagram.privateEntity().child(cur_node, 0);
                j = nbjobs - cur_node.row();
            }
        }

        return sol;
    }
};

struct WeightBDDdouble : WeightBDD<WeightBDDdouble, double> {
    WeightBDDdouble(double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : WeightBDD<WeightBDDdouble, double>(_pi, _interval_list, _nbjobs){};
};

struct WeightZDDdouble : WeightZDD<WeightZDDdouble, double> {
    WeightZDDdouble(double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : WeightZDD<WeightZDDdouble, double>(_pi, _interval_list, _nbjobs){};
};

struct FarkasZDDdouble : FarkasZDD<FarkasZDDdouble, double> {
    FarkasZDDdouble(
        double *_pi, Job *_jobarray, int _nbjobs, int Hmin, int Hmax)
        : FarkasZDD<FarkasZDDdouble, double>(
              _pi, _jobarray, _nbjobs, Hmin, Hmax){};
};
