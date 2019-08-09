#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <memory>
#include <OptimalSolution.hpp>
#include "NodeBase.hpp"
#include "Label.hpp"

template<typename T = double>
class NodeBdd : public NodeBase
{
    private:
        int weight;

    public:
        Label<NodeBdd<T>,T> forward_label[2];
        Label<NodeBdd<T>,T> backward_label[2];

        NodeBdd<T>* child[2];

        bool calc_yes;
        int key;
        int high_edge_key;
        int low_edge_key;

        /**
         * Constructor
         */
        NodeBdd():
            NodeBase(),
            weight(0),
            forward_label{Label<NodeBdd<T>,T>(this), Label<NodeBdd,T>(this)},
            backward_label{Label<NodeBdd<T>,T>(this), Label<NodeBdd,T>(this)},
            calc_yes(true),
            key(-1),
            high_edge_key(-1),
            low_edge_key(-1)
        {
            child[0] = nullptr;
            child[1] = nullptr;
        };

        NodeBdd(int& _weight, int& _num_layer, bool& _root_node, bool& _terminal_node):
            NodeBase(_num_layer, _root_node, _terminal_node),
            weight(_weight),
            forward_label{Label<NodeBdd<T>,T>(this), Label<NodeBdd,T>(this)},
            backward_label{Label<NodeBdd<T>,T>(this), Label<NodeBdd,T>(this)},
            calc_yes(true),
            key(-1),
            high_edge_key(-1),
            low_edge_key(-1)
        {
            child[0] = nullptr;
            child[1] = nullptr;
        }

        void set_head_node()
        {
            forward_label[0].set_head_node(this);
            forward_label[1].set_head_node(this);
            backward_label[0].set_head_node(this);
            backward_label[1].set_head_node(this);
        }

        NodeBdd(int i, int j) :
            NodeBase(i, j),
            weight(0),
            forward_label{Label<NodeBdd<T>,T>(this), Label<NodeBdd,T>(this)},
            backward_label{Label<NodeBdd<T>,T>(this), Label<NodeBdd,T>(this)},
            calc_yes(true),
            key(-1),
            high_edge_key(-1),
            low_edge_key(-1)
        {
            child[0] = nullptr;
            child[1] = nullptr;
        }

        NodeBdd<T>(const NodeBdd<T>& src) = default;
        NodeBdd<T>(NodeBdd<T>&& src) = default;
        NodeBdd<T>& operator=(const NodeBdd<T>& src) = default;
        NodeBdd<T>& operator=(NodeBdd<T>&& src) = default;

        void set_weight(int _weight)
        {
            weight = _weight;
        }

        int get_weight()
        {
            return weight;
        }

        bool operator!=(NodeBdd const& o) const
        {
            return !operator==(o);
        }

        friend std::ostream& operator<<(std::ostream& os, NodeBdd const& o)
        {
            os << "(" << o.branch[0];

            for (int i = 1; i < 2; ++i) {
                os << "," << o.branch[i];
            }

            return os << ")";
        }

        NodeBdd<T>* init_node(int _weight, bool _root_node = false, bool _terminal_node = false)
        {
            if (!_terminal_node) {
                weight = _weight;
                NodeBase::set_root_node(_root_node);
            } else {
                NodeBase::set_job(nullptr);
                weight = -1;
                NodeBase::set_root_node(_root_node);
                NodeBase::set_terminal_node(_terminal_node);
            }

            return this;
        }

        friend bool operator<(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs)
        {
            return lhs.forward_label[0].f < rhs.forward_label[0].f;
        }

        friend bool operator> (const NodeBdd<T>& lhs, const NodeBdd<T>& rhs)
        {
            return rhs < lhs;
        }
        friend bool operator<=(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs)
        {
            return !(lhs > rhs);
        }
        friend bool operator>=(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs)
        {
            return !(lhs < rhs);
        }
};


#endif // NODE_DURATION_HPP
