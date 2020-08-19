#ifndef ZDD_NODE_HPP
#define ZDD_NODE_HPP
#include <OptimalSolution.hpp>
#include <memory>
#include <vector>
#include "Label.hpp"
#include "NodeBase.hpp"

template <typename T = double>
class NodeZdd;

template <typename T = double>
class SubNodeZdd {
   public:
    int weight;

    Label<SubNodeZdd<T>, T> forward_label[2];
    Label<SubNodeZdd<T>, T> backward_label[2];

    std::shared_ptr<SubNodeZdd<T>> y;
    std::shared_ptr<SubNodeZdd<T>> n;

    bool        calc_yes;
    int         key;
    int         high_edge_key;
    int         low_edge_key;
    NodeId      node_id;
    NodeZdd<T>* node_ptr;

    SubNodeZdd()
        : weight{0},
          forward_label{},
          backward_label{},
          y{nullptr},
          n{nullptr},
          calc_yes{true},
          key{-1},
          high_edge_key{-1},
          low_edge_key{-1},
          node_id{0},
          node_ptr{nullptr} {}

    explicit SubNodeZdd(int _weight)
        : weight{_weight},
          forward_label{},
          backward_label{},
          y{nullptr},
          n{nullptr},
          calc_yes{true},
          key{-1},
          high_edge_key{-1},
          low_edge_key{-1},
          node_id{0},
          node_ptr{nullptr} {}

    explicit SubNodeZdd(int _weight, NodeId _node_id)
        : weight{_weight},
          forward_label{},
          backward_label{},
          y{nullptr},
          n{nullptr},
          calc_yes{true},
          key{-1},
          high_edge_key{-1},
          low_edge_key{-1},
          node_id{_node_id},
          node_ptr{nullptr} {
        for (unsigned i = 0; i < 2; ++i) {
            forward_label[i].set_head_node(this);
            backward_label[i].set_head_node(this);
        }
    }

    explicit SubNodeZdd(int _weight, NodeId _node_id, NodeZdd<T>* _node_ptr)
        : weight{_weight},
          forward_label{},
          backward_label{},
          y{nullptr},
          n{nullptr},
          calc_yes{true},
          key{-1},
          high_edge_key{-1},
          low_edge_key{-1},
          node_id{_node_id},
          node_ptr{_node_ptr} {
        for (unsigned i = 0; i < 2; ++i) {
            forward_label[i].set_head_node(this);
            backward_label[i].set_head_node(this);
        }
    }

    SubNodeZdd(const SubNodeZdd& other) = default;
    SubNodeZdd(SubNodeZdd&& other) = default;
    SubNodeZdd<T>& operator=(const SubNodeZdd<T>& other) = default;
    SubNodeZdd<T>& operator=(SubNodeZdd<T>&& other) = default;

    friend bool operator<(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return lhs.forward_label[0].f < rhs.forward_label[0].f;
    }

    friend bool operator>(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return !(lhs > rhs);
    }
    friend bool operator>=(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return !(lhs < rhs);
    }

    int get_weight() { return weight; }

    Job* get_job() { return node_ptr->get_job(); }
};

template <typename T>
bool compare_sub_nodes(const std::shared_ptr<SubNodeZdd<T>>& lhs,
                       const std::shared_ptr<SubNodeZdd<T>>& rhs) {
    return *lhs < *rhs;
}

template <typename T>
class NodeZdd : public NodeBase {
   public:
    NodeZdd<T>*                                 child[2];
    std::shared_ptr<NodeId>                     ptr_node_id;
    std::vector<std::shared_ptr<SubNodeZdd<T>>> list;
    /**
     * Constructor
     */
    NodeZdd() : NodeBase(), list() {
        child[0] = nullptr;
        child[1] = nullptr;
    };

    // NodeZdd(int& _num_layer, bool& _root_node, bool& _terminal_node)
    //     : NodeBase(_num_layer, _root_node, _terminal_node), list() {
    //     child[0] = nullptr;
    //     child[1] = nullptr;
    // };

    void set_head_node() {
        for (auto& it : list) {
            for (unsigned i = 0; i < 2; ++i) {
                it->node_ptr = this;
            }
        }
    }

    NodeZdd(int i, int j) : NodeBase(i, j), child{nullptr, nullptr}, list() {}

    NodeZdd<T>(const NodeZdd<T>& src) = default;
    NodeZdd<T>(NodeZdd<T>&& src) = default;
    NodeZdd<T>& operator=(const NodeZdd<T>& src) = default;
    NodeZdd<T>& operator=(NodeZdd<T>&& src) = default;

    bool operator==(NodeZdd const& o) const {
        for (int i = 0; i < 2; ++i) {
            if (branch[i] != o.branch[i]) {
                return false;
            }
        }

        return true;
    }

    bool operator!=(NodeZdd const& o) const { return !operator==(o); }

    friend std::ostream& operator<<(std::ostream& os, NodeZdd const& o) {
        os << "(" << o.branch[0];

        for (int i = 1; i < 2; ++i) {
            os << "," << o.branch[i];
        }

        return os << ")";
    }

    void add_sub_node(int                   _weight,
                      NodeId                _node_id,
                      [[maybe_unused]] bool _root_node = false,
                      [[maybe_unused]] bool _terminal_node = false) {
        // if(!_terminal_node) {
        list.push_back(std::make_shared<SubNodeZdd<>>(_weight, _node_id, this));
        // printf("test test %p\n", this);
        // set_root_node(_root_node);
        // set_terminal_node(_terminal_node);
        // } else {
        //     list.push_back(std::make_shared<SubNodeZdd<>>(_weight,_node_id));
        //     set_root_node(_root_node);
        //     set_terminal_node(_terminal_node);
        // }
    }

    std::shared_ptr<SubNodeZdd<>> add_weight(int _weight, NodeId _node_id) {
        for (auto& it : list) {
            if (it->weight == _weight) {
                return it;
            }
        }

        add_sub_node(_weight, _node_id);

        return list.back();
    }

    void set_node_id(NodeId _node_id) {
        for (auto& it : list) {
            it->node_id = _node_id;
            it->node_ptr = this;
        }
    }
};

#endif  // ZDD_NODE_HPP
