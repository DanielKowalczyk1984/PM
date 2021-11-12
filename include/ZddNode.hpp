// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef ZDD_NODE_HPP
#define ZDD_NODE_HPP
#include <array>                  // for array
#include <cstddef>                // for size_t
#include <memory>                 // for shared_ptr, make_shared
#include <ostream>                // for operator<<, ostream
#include <vector>                 // for vector
#include "Job.h"                  // for Job
#include "ModernDD/NodeBase.hpp"  // for NodeBase
#include "ModernDD/NodeId.hpp"    // for NodeId
template <typename N>
class Label;
template <typename T = double>
class NodeZdd;

template <typename T = double>
class SubNodeZdd {
   public:
    int weight{};

    std::array<Label<SubNodeZdd<T>>, 2> forward_label{};
    std::array<Label<SubNodeZdd<T>>, 2> backward_label{};

    std::shared_ptr<SubNodeZdd<T>> y{};
    std::shared_ptr<SubNodeZdd<T>> n{};

    bool        calc_yes{true};
    int         key{-1};
    int         high_edge_key{-1};
    int         low_edge_key{-1};
    NodeId      node_id{};
    NodeZdd<T>* node_ptr{nullptr};

    SubNodeZdd() = default;
    ~SubNodeZdd() = default;

    explicit SubNodeZdd(int _weight) : weight{_weight} {}

    explicit SubNodeZdd(int _weight, NodeId _node_id)
        : weight{_weight},
          node_id{_node_id} {
        // for (unsigned i = 0; i < 2; ++i) {
        //     forward_label[i].set_head_node(this);
        //     backward_label[i].set_head_node(this);
        // }
    }

    explicit SubNodeZdd(int _weight, NodeId _node_id, NodeZdd<T>* _node_ptr)
        : weight{_weight},
          node_id{_node_id},
          node_ptr{_node_ptr} {
        // for (unsigned i = 0; i < 2; ++i) {
        //     forward_label[i].set_head_node(this);
        //     backward_label[i].set_head_node(this);
        // }
    }

    SubNodeZdd(const SubNodeZdd& other) = default;
    SubNodeZdd(SubNodeZdd&& other) noexcept = default;
    SubNodeZdd<T>& operator=(const SubNodeZdd<T>& other) = default;
    SubNodeZdd<T>& operator=(SubNodeZdd<T>&& other) noexcept = default;

    friend bool operator<(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return lhs.forward_label[0].get_f() < rhs.forward_label[0].get_f();
    }

    friend bool operator>(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return rhs < lhs;
    }
    friend bool operator<=(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return lhs <= rhs;
    }
    friend bool operator>=(const SubNodeZdd<T>& lhs, const SubNodeZdd<T>& rhs) {
        return lhs >= rhs;
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
    Job* job{nullptr};

   public:
    std::array<NodeZdd<T>*, 2>                  child{};
    std::shared_ptr<NodeId>                     ptr_node_id;
    std::vector<std::shared_ptr<SubNodeZdd<T>>> list{};
    /**
     * Constructor
     */
    NodeZdd() = default;

    NodeZdd(size_t i, size_t j) : NodeBase(i, j) {}

    NodeZdd(const NodeZdd<T>& src) = default;
    NodeZdd(NodeZdd<T>&& src) noexcept = default;
    NodeZdd<T>& operator=(const NodeZdd<T>& src) = default;
    NodeZdd<T>& operator=(NodeZdd<T>&& src) noexcept = default;
    ~NodeZdd() = default;

    bool operator==(NodeZdd const& o) const {
        for (int i = 0; i < 2; ++i) {
            if ((*this)[i] != o[i]) {
                return false;
            }
        }

        return true;
    }

    void set_job(Job* _job) { job = _job; }

    [[nodiscard]] Job* get_job() const { return job; }

    [[nodiscard]] inline int get_nb_job() const { return job->job; }

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

    void set_node_id_label([[maybe_unused]] NodeId _node_id) {}
};

#endif  // ZDD_NODE_HPP
