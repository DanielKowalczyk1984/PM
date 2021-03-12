#ifndef NODE_BASE_HPP
#define NODE_BASE_HPP

// #include <boost/container_hash/hash_fwd.hpp>
#include <array>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include "Job.h"
#include "NodeId.hpp"
class NodeBase : public std::array<NodeId, 2> {
   private:
    Job* job{nullptr};

   public:
    /**
     * Constructor
     */
    NodeBase() : std::array<NodeId, 2>{} {};

    NodeBase(size_t i, size_t j) : std::array<NodeId, 2>{i, j} {}

    NodeBase(NodeId f0, NodeId f1) : std::array<NodeId, 2>{f0, f1} {}

    NodeBase(const NodeBase& src) = default;
    NodeBase(NodeBase&& src) = default;
    NodeBase& operator=(const NodeBase& src) = default;
    NodeBase& operator=(NodeBase&& src) = default;
    ~NodeBase() = default;

    void set_job(Job* _job) {
        job = _job;
        // if (_job != nullptr) {
        //     variable_key[1].set_j(job->job);
        // }
    }

    [[nodiscard]] Job* get_job() const { return job; }

    [[nodiscard]] inline int get_nb_job() const { return job->job; }

    [[nodiscard]] size_t hash() const {
        // size_t h = branch[0].code();
        size_t h = 0;
        // for (int i = 1; i < 2; ++i) {
        // h = h * 314159257 + branch[1].code() * 271828171;
        for (auto const& it : *this) {
            boost::hash_combine(h, it.code());
        }
        // }
        return h;
    }

    bool operator==(NodeBase const& o) const {
        if (*this != o) {
            return false;
        }
        return true;
    }

    bool operator!=(NodeBase const& o) const { return !operator==(o); }

    friend std::ostream& operator<<(std::ostream& os, NodeBase const& o) {
        os << "(" << o[0];
        for (int i = 1; i < 2; ++i) {
            os << "," << o[i];
        }
        return os << ")";
    }
};

struct InitializedNode : NodeBase {
    InitializedNode() : NodeBase(0, 0) {}

    InitializedNode(NodeId f0, NodeId f1) : NodeBase(f0, f1) {}
};

#endif  // NODE_BASE_HPP
