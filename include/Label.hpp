#ifndef LABEL_HPP
#define LABEL_HPP

#include <limits>               // for numeric_limits
#include "ModernDD/NodeId.hpp"  // for NodeId
struct Job;

template <typename N>
class Label {
   private:
    int       weight{-1};
    bool      high{false};
    Label<N>* prev_label{nullptr};
    NodeId    node_id{};

    double f{std::numeric_limits<double>::max()};
    Job*   label_job{nullptr};
    Job*   prev_job{nullptr};

   public:
    /**
     * Constructor
     */
    Label() = default;
    Label(const Label<N>& src) = default;
    Label(Label<N>&& src) noexcept = default;
    Label<N>& operator=(const Label<N>& src) = default;
    Label<N>& operator=(Label<N>&& src) noexcept = default;
    ~Label() = default;

    void set_f(double _f) { f = _f; }

    void set_job(Job* _job) { label_job = _job; };

    void set_node_id(NodeId _id) { node_id = _id; }

    void reset() {
        f = std::numeric_limits<double>::max();
        prev_label = nullptr;
        high = false;
    }

    [[nodiscard]] double get_f() const { return f; }

    double& get_f() { return f; }

    Label<N>* get_previous() { return prev_label; }

    Job* prev_job_backward() { return prev_job; }

    bool get_high() { return high; }

    Job* get_job() { return label_job; }

    int get_weight() { return weight; }

    NodeId& get_node_id() { return node_id; }

    Job* prev_job_forward() {
        return get_previous() == nullptr ? nullptr : get_previous()->get_job();
    }

    void forward_update(double _f, Label<N>* _prev, bool _high = false) {
        f = _f;
        prev_label = _prev;
        high = _high;
    }

    void forward_update(double _f, Label<N>& _node) {
        f = _f;
        prev_label = _node.prev_label;
        high = _node.high;
    }

    void forward_update(Label<N>& _node) {
        f = _node.f;
        prev_label = _node.prev_label;
        high = _node.high;
    }

    void backward_update(double _f, bool _high = false) {
        f = _f;
        high = _high;
    }

    void backward_update(Label<N>* _n, double _f = .0, bool _high = false) {
        f = _f;
        prev_job = _high ? get_job() : _n->prev_job;
        high = _high;
        prev_label = _n;
    }
};

#endif  // LABEL_HPP
