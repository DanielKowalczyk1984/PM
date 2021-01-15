#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <array>
#include <boost/dynamic_bitset.hpp>
#include <memory>
#include "Label.hpp"
#include "NodeBase.hpp"
#include "gurobi_c++.h"

template <typename T = double>
class NodeBdd : public NodeBase {
   private:
    int weight{};

   public:
    std::array<Label<NodeBdd<T>, T>, 2>                 forward_label{};
    std::array<Label<NodeBdd<T>, T>, 2>                 backward_label{};
    std::array<std::vector<std::weak_ptr<BddCoeff>>, 2> coeff_list{};
    std::array<std::vector<std::weak_ptr<NodeId>>, 2>   in_edges{};

    std::shared_ptr<NodeId> ptr_node_id{nullptr};
    std::array<double, 2>   cost{0.0, 0.0};
    std::array<double, 2>   reduced_cost{0.0, 0.0};
    std::array<double, 2>   lp_x{0.0, 0.0};

    bool                    calc_yes{true};
    bool                    calc_no{true};
    int                     key{-1};
    int                     high_edge_key{-1};
    int                     low_edge_key{-1};
    bool                    visited{false};
    bool                    lp_visited{false};
    boost::dynamic_bitset<> all{};
    std::array<double, 2>   backward_distance{};
    std::array<int, 2>      in_degree{};
    std::array<GRBVar, 2>   y;
    std::array<GRBVar, 2>   r;
    GRBVar                  sigma;
    std::array<double, 2>   coeff_cut{0.0, 0.0};

    /**
     * Constructor
     */
    NodeBdd() = default;

    // void set_head_node() {
    //     forward_label[0].set_head_node(this);
    //     forward_label[1].set_head_node(this);
    //     backward_label[0].set_head_node(this);
    //     backward_label[1].set_head_node(this);
    // }

    NodeBdd(int i, int j) : NodeBase(i, j) {}

    NodeBdd<T>(const NodeBdd<T>& src) = default;
    NodeBdd<T>(NodeBdd<T>&& src) noexcept = default;
    NodeBdd<T>& operator=(const NodeBdd<T>& src) = default;
    NodeBdd<T>& operator=(NodeBdd<T>&& src) noexcept = default;

    void set_weight(int _weight) { weight = _weight; }

    [[nodiscard]] int get_weight() const { return weight; }

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

    void init_node(int                   _weight,
                   [[maybe_unused]] bool _root_node = false,
                   bool                  _terminal_node = false) {
        if (!_terminal_node) {
            weight = _weight;
        } else {
            NodeBase::set_job(nullptr);
            weight = -1;
        }
    }

    void set_node_id_label(NodeId _id) {
        for (int j = 0; j < 2; j++) {
            backward_label[j].set_node_id(_id);
            forward_label[j].set_node_id(_id);
        }
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
