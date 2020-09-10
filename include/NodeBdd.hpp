#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <boost/dynamic_bitset.hpp>
#include <memory>
#include "Label.hpp"
#include "NodeBase.hpp"
#include "gurobi_c++.h"

template <typename T = double>
class NodeBdd : public NodeBase {
   private:
    int weight;

   public:
    Label<NodeBdd<T>, T>                 forward_label[2];
    Label<NodeBdd<T>, T>                 backward_label[2];
    std::vector<std::weak_ptr<BddCoeff>> coeff_list[2];
    std::vector<std::weak_ptr<NodeId>>   in_edges[2];

    NodeBdd<T>*             child[2];
    std::shared_ptr<NodeId> ptr_node_id;
    double                  cost[2];
    double                  reduced_cost[2];
    double                  lp_x[2];

    bool                    calc_yes;
    bool                    calc_no;
    int                     key;
    int                     high_edge_key;
    int                     low_edge_key;
    bool                    visited;
    bool                    lp_visited;
    int                     lp_key;
    boost::dynamic_bitset<> all;
    int                     backward_distance[2];
    int                     in_degree_0;
    int                     in_degree_1;
    GRBVar                  y[2];
    GRBVar                  r[2];
    double                  coeff_cut[2];

    /**
     * Constructor
     */
    NodeBdd()
        : NodeBase(),
          weight(0),
          forward_label{Label<NodeBdd<T>, T>(this), Label<NodeBdd, T>(this)},
          backward_label{Label<NodeBdd<T>, T>(this), Label<NodeBdd, T>(this)},
          ptr_node_id(nullptr),
          cost{0.0, 0.0},
          reduced_cost{0.0, 0.0},
          lp_x{0.0, 0.0},
          calc_yes(true),
          calc_no(true),
          key(-1),
          high_edge_key(-1),
          low_edge_key(-1),
          visited(false),
          lp_visited(false),
          lp_key(-1),
          in_degree_0(0),
          in_degree_1(0),
          coeff_cut{0.0, 0.0} {
        child[0] = nullptr;
        child[1] = nullptr;
    };

    NodeBdd(int&  _weight,
            int&  _num_layer,
            bool& _root_node,
            bool& _terminal_node)
        : NodeBase(_num_layer, _root_node, _terminal_node),
          weight(_weight),
          forward_label{Label<NodeBdd<T>, T>(this), Label<NodeBdd, T>(this)},
          backward_label{Label<NodeBdd<T>, T>(this), Label<NodeBdd, T>(this)},
          ptr_node_id(nullptr),
          cost{0.0, 0.0},
          reduced_cost{0.0, 0.0},
          lp_x{0.0, 0.0},
          calc_yes(true),
          calc_no(true),
          key(-1),
          high_edge_key(-1),
          low_edge_key(-1),
          visited(false),
          lp_visited(false),
          lp_key(-1),
          in_degree_0(0),
          in_degree_1(0),
          coeff_cut{0.0, 0.0} {
        child[0] = nullptr;
        child[1] = nullptr;
    }

    void set_head_node() {
        forward_label[0].set_head_node(this);
        forward_label[1].set_head_node(this);
        backward_label[0].set_head_node(this);
        backward_label[1].set_head_node(this);
    }

    NodeBdd(int i, int j)
        : NodeBase(i, j),
          weight(0),
          forward_label{Label<NodeBdd<T>, T>(this), Label<NodeBdd, T>(this)},
          backward_label{Label<NodeBdd<T>, T>(this), Label<NodeBdd, T>(this)},
          ptr_node_id(nullptr),
          cost{0.0, 0.0},
          reduced_cost{0.0, 0.0},
          lp_x{0.0, 0.0},
          calc_yes(true),
          calc_no(true),
          key(-1),
          high_edge_key(-1),
          low_edge_key(-1),
          visited(false),
          lp_visited(false),
          lp_key(-1),
          in_degree_0(0),
          in_degree_1(0),
          coeff_cut{0.0, 0.0} {
        child[0] = nullptr;
        child[1] = nullptr;
    }

    NodeBdd<T>(const NodeBdd<T>& src) = default;
    NodeBdd<T>(NodeBdd<T>&& src) = default;
    NodeBdd<T>& operator=(const NodeBdd<T>& src) = default;
    NodeBdd<T>& operator=(NodeBdd<T>&& src) = default;

    void set_weight(int _weight) { weight = _weight; }

    int get_weight() const { return weight; }

    void reset_reduced_costs() {
        reduced_cost[0] = 0.0;
        reduced_cost[1] = cost[1];
    }

    void reset_reduced_costs_farkas() {
        reduced_cost[0] = reduced_cost[1] = 0.0;
    }

    void reset_lp_x() {
        lp_x[0] = lp_x[1] = 0.0;
        lp_visited = false;
        lp_key = -1;
    }

    void add_coeff_list(std::shared_ptr<BddCoeff>& ptr, int high) {
        coeff_list[high].push_back(ptr);
    }

    void adjust_reduced_costs(double _x, int i) { reduced_cost[i] -= _x; }

    bool operator!=(NodeBdd const& o) const { return !operator==(o); }

    friend std::ostream& operator<<(std::ostream& os, NodeBdd const& o) {
        os << "(" << o.branch[0];

        for (int i = 1; i < 2; ++i) {
            os << "," << o.branch[i];
        }

        return os << ")";
    }

    NodeBdd<T>* init_node(int                   _weight,
                          [[maybe_unused]] bool _root_node = false,
                          bool                  _terminal_node = false) {
        if (!_terminal_node) {
            weight = _weight;
        } else {
            NodeBase::set_job(nullptr);
            weight = -1;
        }

        return this;
    }

    friend bool operator<(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs) {
        return lhs.forward_label[0].f < rhs.forward_label[0].f;
    }

    friend bool operator>(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs) {
        return !(lhs > rhs);
    }
    friend bool operator>=(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs) {
        return !(lhs < rhs);
    }
};

#endif  // NODE_DURATION_HPP
