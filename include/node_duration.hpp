#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP
#include <memory>
#include <OptimalSolution.hpp>
#include "NodeBase.hpp"
#include "Label.hpp"

template<typename T = double>
class Node : public NodeBase
{
    private:
        int weight;

    public:
        Label<Node<T>,T> forward_label[2];
        Label<Node<T>,T> backward_label[2];

        Node<T>* child[2];

        bool calc_yes;
        int key;
        int high_edge_key;
        int low_edge_key;

        /**
         * Constructor
         */
        Node():
            NodeBase(),
            weight(0),
            forward_label{Label<Node<T>,T>(this), Label<Node,T>(this)},
            backward_label{Label<Node<T>,T>(this), Label<Node,T>(this)},
            calc_yes(true),
            key(-1),
            high_edge_key(-1),
            low_edge_key(-1)
        {
            child[0] = nullptr;
            child[1] = nullptr;
        };

        Node(int& _weight, int& _num_layer, bool& _root_node, bool& _terminal_node):
            NodeBase(_num_layer, _root_node, _terminal_node),
            weight(_weight),
            forward_label{Label<Node<T>,T>(this), Label<Node,T>(this)},
            backward_label{Label<Node<T>,T>(this), Label<Node,T>(this)},
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

        Node(int i, int j) :
            NodeBase(i, j),
            weight(0),
            forward_label{Label<Node<T>,T>(this), Label<Node,T>(this)},
            backward_label{Label<Node<T>,T>(this), Label<Node,T>(this)},
            calc_yes(true),
            key(-1),
            high_edge_key(-1),
            low_edge_key(-1)
        {
            child[0] = nullptr;
            child[1] = nullptr;
        }

        Node<T>(const Node<T>& src) = default;
        Node<T>(Node<T>&& src) = default;
        Node<T>& operator=(const Node<T>& src) = default;
        Node<T>& operator=(Node<T>&& src) = default;

        void set_weight(int _weight)
        {
            weight = _weight;
        }

        int get_weight()
        {
            return weight;
        }

        bool operator!=(Node const& o) const
        {
            return !operator==(o);
        }

        friend std::ostream& operator<<(std::ostream& os, Node const& o)
        {
            os << "(" << o.branch[0];

            for (int i = 1; i < 2; ++i) {
                os << "," << o.branch[i];
            }

            return os << ")";
        }

        Node<T>* init_node(int _weight, bool _root_node = false, bool _terminal_node = false)
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

        friend bool operator<(const Node<T>& lhs, const Node<T>& rhs)
        {
            return lhs.forward_label[0].f < rhs.forward_label[0].f;
        }

        friend bool operator> (const Node<T>& lhs, const Node<T>& rhs)
        {
            return rhs < lhs;
        }
        friend bool operator<=(const Node<T>& lhs, const Node<T>& rhs)
        {
            return !(lhs > rhs);
        }
        friend bool operator>=(const Node<T>& lhs, const Node<T>& rhs)
        {
            return !(lhs < rhs);
        }
};


#endif // NODE_DURATION_HPP
