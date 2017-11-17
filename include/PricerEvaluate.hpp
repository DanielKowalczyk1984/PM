#include <solution.h>
#include <interval.h>
#include <boost/dynamic_bitset.hpp>
#include <cfloat>
#include <climits>
#include <tdzdd/DdEval.hpp>
#include <tdzdd/dd/NodeTable.hpp>
#include <gurobi_c++.h>

void reset_total();
template<typename T>
class edge;

using namespace std;
template <typename T>
class node {
public:
    int      layer;
    int      weight;
    T obj, obj0, obj1;
    T b;
    T c;
    T dist;
    bool calc;
    bool calc0;
    bool remove_node;
    bool take;
    Job *prev_job;
    shared_ptr<node<T>> prev_node;
    bool terminal;
    shared_ptr<node<T>> y;
    shared_ptr<node<T>> n;
    vector<weak_ptr<edge<T>>> out_edge;
    vector<weak_ptr<edge<T>>> in_edge;

    /** Default Constructor */
    node()
        : layer(0),
          weight(0),
          obj(0.0),
          obj0(0.0),
          obj1(0.0),
          b(0.0),
          c(0.0),
          dist(-DBL_MAX),
          calc(true),
          calc0(true),
          remove_node(false),
          take(false),
          prev_job(nullptr),
          prev_node(nullptr),
          terminal(false),
          y(nullptr),
          n(nullptr)
    {
        // out_edge.reserve(2);
        // in_edge.reserve(2);
    };

    /** copy constructor */
    node<T>(const node<T>& other) :
          layer(other.layer),
          weight(other.weight),
          obj(other.obj),
          obj0(other.obj0),
          obj1(other.obj1),
          b(other.b),
          c(other.c),
          dist(other.dist),
          calc(other.calc),
          calc0(other.calc0),
          remove_node(other.remove_node),
          take(other.take),
          prev_job(other.prev_job),
          prev_node(other.prev_job),
          terminal(other.terminal),
          y(other.y),
          n(other.n),
          out_edge(other.out_edge),
          in_edge(other.in_edge) {

    }

    /** move constructor */
    node<T>(node<T>&& other) noexcept:
    layer(other.layer),
    weight(other.weight),
    obj(other.obj),
    obj0(other.obj0),
    obj1(other.obj1),
    b(other.b),
    c(other.c),
    dist(other.dist),
    calc(other.calc),
    calc0(other.calc0),
    remove_node(other.remove_node),
    take(other.take),
    prev_job(other.prev_job),
    prev_node(other.prev_job),
    terminal(other.terminal),
    y(other.y),
    n(other.n),
    out_edge(other.out_edge),
    in_edge(other.in_edge){
        other.y.reset();
        other.n.reset();
        other.prev_node.reset();
        other.prev_job = nullptr;
        other.in_edge.clear();
        other.out_edge.clear();
    };

    /** copy assignment */
    node<T>& operator=(const node<T> other){
        node<T> tmp(other);
        *this = move(tmp);
        return *this;
    }

    /** move assignment */
    node<T>& operator= (node<T>&& other) noexcept{
        if(this != &other) {
            layer =other.layer;
                weight =other.weight;
                obj =other.obj;
                obj0 =other.obj0;
                obj1 =other.obj1;
                b =other.b;
                c =other.c;
                dist =other.dist;
                calc =other.calc;
                calc0 =other.calc0;
                remove_node =other.remove_node;
                take =other.take;
                prev_job =other.prev_job;
                prev_node =other.prev_job;
                terminal =other.terminal;
                y =other.y;
                n =other.n;
                other.y = nullptr;
                other.n = nullptr;
                other.prev_node = nullptr;
                other.prev_job = nullptr;
                other.in_edge.clear();
                other.out_edge.clear();
        }
        return *this;
    }

    ~node() noexcept {}  ;
};


template<typename T>
class edge {
private:
    static int _total;


public:
    int id;
    T cost;
    Job *job;
    shared_ptr<node<T>> out;
    shared_ptr<node<T>> in;
    GRBVar v;



    /** Default constructor */
    edge() : id(0), cost(0.0), job (nullptr), out(nullptr), in(nullptr)
    {
    }

    edge(T _cost, Job *_job, shared_ptr<node<T>>& _out, shared_ptr<node<T>>& _in) : cost(_cost), job(_job), out(_out), in(_in)
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
        *this = move(tmp); // re-use move-assignment
        return *this;
    }

    /** Move assignment operator */
    edge& operator= (edge&& other) noexcept
    {
        if(this != &other) {
            id = other.id;
            cost = other.cost;
            out = other.out;
            in = other.in;
            job = other.job;
            v = other.v;
            other.in = nullptr;
            other.out = nullptr;
            other.job = nullptr;
        }
        return *this;
    }

    /** Destructor */
    ~edge() noexcept{
    }

    friend ostream& operator<<(ostream&              os,
                                    edge<T> const& o)
    {
        os << "id = " << o.id << "\n";
        return os;
    };

    friend void reset_total(){
        _total = 0;
    }

};

template<typename T>
int edge<T>::_total = 0;

template <typename T>
class PricerWeightZDD {
public:
    vector<shared_ptr<node<T>>> list;

    PricerWeightZDD() {};

    ~PricerWeightZDD() noexcept{
        list.clear();
    }

    shared_ptr<node<T>>& add_weight(int _weight, int layer, int njobs)
    {
        for (auto& it : list) {
            if (it->weight == _weight || it->terminal) {
                return it;
            }
        }

        shared_ptr<node<T>> p = make_shared<node<T>>();
        p->layer = layer;
        p->weight = _weight;
        list.push_back(p);
        return (list.back());
    }

    void add_terminal_node(int j, int layer, int njobs){
        shared_ptr<node<T>> p = make_shared<node<T>>();
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

/**
 * BDD
 */
template<typename T>
class PricerInfoBDD {
  public:
    T obj;
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

template <typename T>
class Optimal_Solution {
public:
    T                obj;
    int              cost;
    int              C_max;
    GPtrArray       *jobs;
    GPtrArray       *e_list;


    /** Default constructor */
    Optimal_Solution() : obj(0), cost(0), C_max(0), jobs(g_ptr_array_new()), e_list(g_ptr_array_new())
    {
    }

    /** Copy constructor */
    Optimal_Solution(const Optimal_Solution& other) :
        obj(other.obj), cost(other.cost), C_max(other.C_max),
        jobs(g_ptr_array_sized_new(other.jobs->len)),
        e_list(g_ptr_array_sized_new(other.e_list->len))
    {
        for (unsigned i = 0; i < other.jobs->len; ++i) {
            g_ptr_array_add(jobs, g_ptr_array_index(other.jobs, i));
        }

        for (unsigned i = 0; i < other.edges->len; ++i) {
            g_ptr_array_add(jobs, g_ptr_array_index(other.edges, i));
        }
    }

    /** Move constructor */
    Optimal_Solution(Optimal_Solution&& other)
    noexcept :  /* noexcept needed to enable optimizations in containers */
        obj(other.obj), cost(other.cost), C_max(other.C_max), jobs(other.jobs), e_list(other.e_list)
    {
        other.jobs = nullptr;
        other.e_list = nullptr;
    }

    /** Copy assignment operator */
    Optimal_Solution& operator= (const Optimal_Solution& other)
    {
        Optimal_Solution tmp(other);         // re-use copy-constructor
        *this = move(tmp); // re-use move-assignment
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
        g_ptr_array_free(e_list, TRUE);
        e_list = other.e_list;
        other.e_list = nullptr;
        return *this;
    }

    Optimal_Solution(PricerInfoBDD<T> *node) {
        cost = 0;
        obj = node->obj;
        C_max = 0;
        jobs = g_ptr_array_sized_new(node->jobs.size());
        for(const auto& it: node->jobs){
            g_ptr_array_add(jobs, it);
            C_max += ((Job *) it)->processingime;
            cost += value_Fj(C_max, it);
        }
    }

    /** Destructor */
    ~Optimal_Solution() noexcept
    {
        if (jobs) {
            g_ptr_array_free(jobs, TRUE);
        }

        if(e_list) {
            g_ptr_array_free(e_list, TRUE);
        }
    }

    friend ostream& operator<<(ostream&              os,
                                    Optimal_Solution<T> const& o)
    {
        os << "obj = " << o.obj << "," << endl
           << "cost = " << o.cost << " C_max = " << o.C_max << endl;
        g_ptr_array_foreach(o.jobs, g_print_machine, NULL);
        return os;
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

    void evalNode(PricerInfoBDD<T> &n, int i,
                  tdzdd::DdValues<PricerInfoBDD<T>, 2>    &values) const {
        int j = nlayers - i;
        assert(j >= 0 && j <= nlayers - 1);
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          interval_list, j);
        Job *tmp_j = tmp_pair->j;
        PricerInfoBDD<T> *n0 = values.get_ptr(0);
        PricerInfoBDD<T> *n1 = values.get_ptr(1);


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

    Optimal_Solution<T> get_objective(PricerInfoBDD<T> *n) {
        Optimal_Solution<T> sol(n);
        return sol;
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
            shared_ptr<node<T>> p0 = it->n;
            shared_ptr<node<T>> p1 = it->y;
            shared_ptr<node<T>> aux;
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
        shared_ptr<node<T>> ptr = (*data_table)[f->row()][f->col()].list[sol.C_max];
        job_interval_pair *tmp_pair;
        Job *tmp_j;


        while (ptr->layer != nlayers) {
            tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, ptr->layer);
            tmp_j = tmp_pair->j;

            if (ptr->obj == ptr->obj1) {
                g_ptr_array_add(sol.jobs, tmp_j);
                auto e = ptr->out_edge[0].lock();
                g_ptr_array_add(sol.e_list, &(e->id));
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                ptr = ptr->y;
            } else {
                auto e = ptr->out_edge[1].lock();
                g_ptr_array_add(sol.e_list, &(e->id));
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

struct DurationBDDdouble : DurationBDD<DurationBDDdouble, double> {
    DurationBDDdouble(double *_pi, GPtrArray *_interval_list, int _nbjobs)
        : DurationBDD<DurationBDDdouble, double>(_pi, _interval_list, _nbjobs) {};
};

