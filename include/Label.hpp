#ifndef LABEL_HPP
#define LABEL_HPP

#include <limits>      // for numeric_limits
#include "ModernDD/NodeId.hpp"  // for NodeId
struct Job;

template <typename N, typename T>
class Label {
   private:
    int          weight{-1};
    bool         high{false};
    Label<N, T>* prev_label{nullptr};
    NodeId       node_id{};

    T    f{std::numeric_limits<double>::max()};
    Job* label_job{nullptr};
    Job* prev_job{nullptr};

   public:
    /**
     * Constructor
     */
    Label() = default;
    Label(const Label<N, T>& src) = default;
    Label(Label<N, T>&& src) noexcept = default;
    Label<N, T>& operator=(const Label<N, T>& src) = default;
    Label<N, T>& operator=(Label<N, T>&& src) noexcept = default;
    ~Label() = default;

    void set_f(T _f) { f = _f; }

    void set_job(Job* _job) { label_job = _job; };

    void set_node_id(NodeId _id) { node_id = _id; }

    void reset() {
        f = std::numeric_limits<double>::max();
        prev_label = nullptr;
        high = false;
    }

    T get_f() const { return f; }

    T& get_f() { return f; }

    Label<N, T>* get_previous() { return prev_label; }

    Job* prev_job_backward() { return prev_job; }

    bool get_high() { return high; }

    Job* get_job() { return label_job; }

    int get_weight() { return weight; }

    NodeId& get_node_id() { return node_id; }

    Job* prev_job_forward() {
        return get_previous() == nullptr ? nullptr : get_previous()->get_job();
    }

    void forward_update(T _f, Label<N, T>* _prev, bool _high = false) {
        f = _f;
        prev_label = _prev;
        high = _high;
    }

    void forward_update(T _f, Label<N, T>& _node) {
        f = _f;
        prev_label = _node.prev_label;
        high = _node.high;
    }

    void forward_update(Label<N, T>& _node) {
        f = _node.f;
        prev_label = _node.prev_label;
        high = _node.high;
    }

    void backward_update(T _f, bool _high = false) {
        f = _f;
        high = _high;
    }

    void backward_update(Label<N, T>* _n, T _f = 0, bool _high = false) {
        f = _f;
        prev_job = _high ? get_job() : _n->prev_job;
        high = _high;
        prev_label = _n;
    }
};

#endif  // LABEL_HPP
