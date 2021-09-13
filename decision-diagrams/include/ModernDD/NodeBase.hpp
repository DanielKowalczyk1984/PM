#ifndef NODE_BASE_HPP
#define NODE_BASE_HPP

#include <array>                                // for array
#include <boost/container_hash/extensions.hpp>  // for hash_combine
#include <cstddef>                              // for size_t
#include <ostream>                              // for operator<<, ostream
#include "NodeId.hpp"                           // for NodeId, operator<<

class NodeBase : public std::array<NodeId, 2> {
   public:
    NodeBase() : std::array<NodeId, 2>{} {};
    NodeBase(size_t i, size_t j) : std::array<NodeId, 2>{i, j} {}
    NodeBase(NodeId f0, NodeId f1) : std::array<NodeId, 2>{f0, f1} {}

    NodeBase(const NodeBase& src) = default;
    NodeBase(NodeBase&& src) = default;
    NodeBase& operator=(const NodeBase& src) = default;
    NodeBase& operator=(NodeBase&& src) = default;
    ~NodeBase() = default;

    [[nodiscard]] size_t hash() const {
        size_t h = 0;
        for (auto const& it : *this) {
            boost::hash_combine(h, it.code());
        }
        return h;
    }

    // bool operator==(NodeBase const& o) const { return *this == o; }

    // bool operator!=(NodeBase const& o) const { return !operator==(o); }

    friend std::ostream& operator<<(std::ostream& os, NodeBase const& o) {
        os << "(" << o[0];
        for (auto i = 1UL; i < 2; ++i) {
            os << "," << o[i];
        }
        return os << ")";
    }
};

// struct InitializedNode : NodeBase {
//     InitializedNode() : NodeBase(0, 0) {}
//     InitializedNode(NodeId f0, NodeId f1) : NodeBase(f0, f1) {}
// };

#endif  // NODE_BASE_HPP
