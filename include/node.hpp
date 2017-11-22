#ifndef NODE_HPP
#define NODE_HPP
#include <solution.h>
#include <memory>
#include <vector>
#include <gurobi_c++.h>

void reset_total();
template<typename T>
class edge;

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
    std::shared_ptr<node<T>> prev_node;
    bool terminal;
    std::shared_ptr<node<T>> y;
    std::shared_ptr<node<T>> n;
    std::vector<std::weak_ptr<edge<T>>> out_edge;
    std::vector<std::weak_ptr<edge<T>>> in_edge;

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
    std::shared_ptr<node<T>> out;
    std::shared_ptr<node<T>> in;
    GRBVar v;



    /** Default constructor */
    edge() : id(0), cost(0.0), job (nullptr), out(nullptr), in(nullptr)
    {
    }

    edge(T _cost, Job *_job, std::shared_ptr<node<T>>& _out, std::shared_ptr<node<T>>& _in) : cost(_cost), job(_job), out(_out), in(_in)
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

    friend std::ostream& operator<<(std::ostream&              os,
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


#endif // NODE_HPP


