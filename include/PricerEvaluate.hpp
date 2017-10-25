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
#include <gurobi_c++.h>


template <typename T>
class node {
public:
    int id;
    int      layer;
    int      weight;
    T        obj;
    T b;
    T c;
    T dist;
    bool calc;
    bool calc0;
    bool remove_node;
    bool take;
    Job *prev;
    node<T> *y;
    node<T> *n;
    boost::dynamic_bitset<> x;

    node()
        : layer(0),
          weight(0),
          obj(0.0),
          b(0.0),
          c(0.0),
          calc(true),
          calc0(true),
          remove_node(false),
          take(false),
          prev(nullptr),
          y(nullptr),
          n(nullptr)
    {
    };


    node<T>(const node<T>& other)
    {
        layer = other.layer;
        weight = other.weight;
        obj = other.obj;
        b = other.obj;
        take = other.take;
        y = other.y;
        n = other.n;
    };

    ~node(){};
};

template<typename T>
class edge {
private:
    static int _total;


public:
    int id;
    T cost;
    Job *job;
    node<T> *out;
    node<T> *in;
    GRBVar v;



    /** Default constructor */
    edge() : id(0), cost(0.0), job (nullptr), out(nullptr), in(nullptr)
    {
    }

    edge(T _cost, Job *_job, node<T> *_out, node<T> *_in) : cost(_cost), job(_job), out(_out), in(_in)
    {
        id = _total++;
    }

    /** Copy constructor */
    edge(const edge& other) :
        id(other.id), cost(other.cost), job(other.job), out(other.out), in(other.in), v(other.v)
    {
    }

    /** Move constructor */
    edge(edge&& other)
    noexcept :  /* noexcept needed to enable optimizations in containers */
        id(other.id), cost(other.cost), job(other.job), out(other.out), in(other.in), v(other.v)
    {
        other.job = nullptr;
        other.out = nullptr;
        other.in = nullptr;
    }

    /** Copy assignment operator */
    edge& operator= (const edge& other)
    {
        edge tmp(other);         // re-use copy-constructor
        *this = std::move(tmp); // re-use move-assignment
        return *this;
    }

    /** Move assignment operator */
    edge& operator= (edge&& other) noexcept
    {
        id = other.id;
        cost = other.cost;
        out = other.out;
        in = other.in;
        job = other.job;
        v = other.v;
        other.in = nullptr;
        other.out = nullptr;
        other.job = nullptr;
        return *this;
    }

    /** Destructor */
    ~edge() noexcept
    {
    }

};

template<typename T>
int edge<T>::_total = 0;

template <typename T>
using my_iterator = typename std::vector<node<T> *>::iterator;

template <typename T>
class PricerWeightZDD {
public:
    std::vector<node<T> *> list;

    PricerWeightZDD() {};

    ~PricerWeightZDD()
    {
        for (my_iterator<T> it = list.begin(); it != list.end(); it++) {
            delete *it;
        }
    }

    node<T> *add_weight(int _weight, int layer, int njobs)
    {
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
        p->prev = nullptr;
        p->dist = -DBL_MAX;
        p->x.resize(njobs);
        list.push_back(p);
        return (list.back());
    }

    void init_terminal_node(int j, int njobs)
    {
        for (my_iterator<T> it = list.begin(); it != list.end(); it++) {
            (*it)->obj = j ? 0.0 : -DBL_MAX / 100.0;
            (*it)->take = false;
            (*it)->prev = nullptr;
            (*it)->x.resize(njobs);
        }
    }
};

template <typename T>
class PricerFarkasZDD {
public:
    T    obj;
    bool take;

    PricerFarkasZDD() : obj(0), take(0) {};

    ~PricerFarkasZDD() {};

    void init_terminal_node(int one) { obj = one ? 0.0 : 1871286761.0; }

    void init_node() { take = false; }
};

template <typename T>
class Optimal_Solution {
public:
    T                obj;
    int              cost;
    int              C_max;
    GPtrArray       *jobs;


    /** Default constructor */
    Optimal_Solution() : obj(0), cost(0), C_max(0), jobs(g_ptr_array_new())
    {
    }

    /** Copy constructor */
    Optimal_Solution(const Optimal_Solution& other) :
        obj(other.obj), cost(other.cost), C_max(other.C_max),
        jobs(g_ptr_array_sized_new(other.jobs->len))
    {
        for (unsigned i = 0; i < other.jobs->len; ++i) {
            g_ptr_array_add(jobs, g_ptr_array_index(other.jobs, i));
        }
    }

    /** Move constructor */
    Optimal_Solution(Optimal_Solution&& other)
    noexcept :  /* noexcept needed to enable optimizations in containers */
        obj(other.obj), cost(other.cost), C_max(other.C_max), jobs(other.jobs)
    {
        other.jobs = nullptr;
    }

    /** Copy assignment operator */
    Optimal_Solution& operator= (const Optimal_Solution& other)
    {
        Optimal_Solution tmp(other);         // re-use copy-constructor
        *this = std::move(tmp); // re-use move-assignment
        return *this;
    }

    /** Move assignment operator */
    Optimal_Solution& operator= (Optimal_Solution&& other) noexcept
    {
        obj = other.obj;
        cost = other.cost;
        C_max = other.C_max;
        g_ptr_array_free(jobs, TRUE);
        jobs = other.jobs;
        other.jobs = nullptr;
        return *this;
    }

    /** Destructor */
    ~Optimal_Solution() noexcept
    {
        if (jobs) {
            g_ptr_array_free(jobs, TRUE);
        }
    }

    friend std::ostream& operator<<(std::ostream&              os,
                                    Optimal_Solution<T> const& o)
    {
        os << "obj = " << o.obj << "," << std::endl
           << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;
        g_ptr_array_foreach(o.jobs, g_print_machine, NULL);
        return os;
    };

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
        for (my_iterator<T> it = n.list.begin(); it != n.list.end(); it++) {
            (*it)->obj = pi[nbjobs];
            (*it)->prev = nullptr;
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

        for (my_iterator<T> it = n->list.begin(); it != n->list.end(); it++) {
            int      weight = ((*it)->weight);
            node<T> *p0 = (*it)->n;
            node<T> *p1 = (*it)->y;
            T        obj0 = p0->obj;
            T        obj1 = p1->obj;
            result = - value_Fj(weight + tmp_j->processingime, tmp_j) + pi[tmp_j->job];

            if (((*it)->calc)) {
                if (obj0 <= obj1 + result) {
                    if(1 != p1->x[tmp_j->job]) {
                        (*it)->obj = obj1 + result;
                        (*it)->take = true;
                        (*it)->prev = tmp_j;
                        (*it)->x = p1->x;
                        (*it)->x[tmp_j->job] = 1;
                    } else {
                        (*it)->obj = obj0;
                        (*it)->take = false;
                        (*it)->prev = p0->prev;
                        (*it)->x = p0->x;
                    }
                } else {
                    (*it)->obj = obj0;
                    (*it)->take = false;
                    (*it)->prev = p0->prev;
                    (*it)->x = p0->x;
                }

                if ((*it)->b < obj1 + result) {
                    (*it)->b = obj1 + result;
                }

                if((*it)->c < obj0) {
                    (*it)->c = obj0;
                }
            } else {
                (*it)->obj = obj0;
                (*it)->take = false;
                (*it)->prev = p0->prev;
                (*it)->x = p0->x;
            }
        }
    }

    void initializenode(PricerWeightZDD<T>& n)
    {
        for (my_iterator<T> it = n.list.begin(); it != n.list.end(); it++) {
            (*it)->obj = 0.0;
            (*it)->b = pi[nbjobs];
            (*it)->c = pi[nbjobs];
            (*it)->dist = -DBL_MAX;
            (*it)->take = false;
        }
    }

    Optimal_Solution<T> get_objective(tdzdd::NodeTableHandler<2> diagram,
                                      tdzdd::DataTable<PricerWeightZDD<T>> *data_table,
                                      const tdzdd::NodeId *f)
    {
        Optimal_Solution<T> sol;
        sol.obj = (*data_table)[f->row()][f->col()].list[sol.C_max]->obj;
        node<T> *ptr = (*data_table)[f->row()][f->col()].list[sol.C_max];
        job_interval_pair *tmp_pair;
        Job *tmp_j;

        while (ptr->layer != nlayers) {
            tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, ptr->layer);
            tmp_j = tmp_pair->j;

            if (ptr->take) {
                g_ptr_array_add(sol.jobs, tmp_j);
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                ptr = ptr->y;
            } else {
                ptr = ptr->n;
            }
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
    T   *pi;
    GPtrArray *interval_list;
    int nlayers;
    int nbjobs;
    job_interval_pair *tmp_pair;
    Job *tmp_j;
    interval *tmp_interval;

public:
    FarkasZDD(T *_pi, GPtrArray *_interval_list, int _nbjobs)
        : pi(_pi),
          interval_list(_interval_list),
          nbjobs(_nbjobs) {};

    void evalTerminal(PricerFarkasZDD<T>& n) { n.obj = pi[nbjobs]; }

    void evalNode(PricerFarkasZDD<T> *n,
                  int                 i,
                  tdzdd::DdValues<PricerFarkasZDD<T>, 2>& values) const
    {
        int j = nlayers - i;
        assert(j >= 0 && j <= nbjobs - 1);
        tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
        tmp_j = tmp_pair->j;
        PricerFarkasZDD<T> *n0 = values.get_ptr(0);
        PricerFarkasZDD<T> *n1 = values.get_ptr(1);

        if (n0->obj < n1->obj + pi[tmp_j->job]) {
            n->obj = n0->obj;
            n->take = false;
        } else {
            n->obj = n1->obj + pi[tmp_j->job];
            n->take = true;
        }
    }

    void initializenode(PricerFarkasZDD<T>& n) { n.take = false; }

    Optimal_Solution<T> get_objective(
        tdzdd::NodeTableHandler<2>            diagram,
        tdzdd::DataTable<PricerFarkasZDD<T>> *data_table,
        const tdzdd::NodeId                  *f)
    {
        Optimal_Solution<T> sol;
        sol.obj = (*data_table)[f->row()][f->col()].obj;
        sol.jobs = g_ptr_array_new();
        tdzdd::NodeId cur_node = *f;
        int           j = nlayers - cur_node.row();

        while (cur_node.row() != 0) {
            if ((*data_table)[cur_node.row()][cur_node.col()].take) {
                g_ptr_array_add(sol.jobs, tmp_j);
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                cur_node = diagram.privateEntity().child(cur_node, 1);
            } else {
                cur_node = diagram.privateEntity().child(cur_node, 0);
            }

            j = nlayers - cur_node.row();
        }

        return sol;
    }
};

struct WeightZDDdouble : WeightZDD<WeightZDDdouble, double> {
    WeightZDDdouble(double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : WeightZDD<WeightZDDdouble, double>(_pi, _interval_list, _nbjobs) {};
};

struct FarkasZDDdouble : FarkasZDD<FarkasZDDdouble, double> {
    FarkasZDDdouble(
        double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : FarkasZDD<FarkasZDDdouble, double>(_pi, _interval_list, _nbjobs) {};
};
