#ifndef LABEL_HPP
#define LABEL_HPP

#include <limits>
#include "NodeId.hpp"
#include "job.h"

template <typename N, typename T>
class Label {
   private:
    Label<N, T>* prev{nullptr};
    bool         high{false};
    Job*         job{nullptr};
    int          weight{-1};
    N*           head_node{nullptr};
    NodeId       node_id{};

    T    f{std::numeric_limits<double>::max()};
    Job* prev_job{nullptr};

   public:
    /**
     * Constructor
     */
    Label(T& _f, Label<N, T>*& _prev, bool _high)
        : prev(_prev),
          high(_high),
          f(_f){};

    Label() = default;

    Label<N, T>(const Label<N, T>& src) = default;
    Label<N, T>(Label<N, T>&& src) noexcept = default;
    Label<N, T>& operator=(const Label<N, T>& src) = default;
    Label<N, T>& operator=(Label<N, T>&& src) noexcept = default;
    ~Label<N, T>() = default;

    void set_previous(Label<N, T>* _prev) { prev = _prev; }

    void set_f(T _f) { f = _f; }

    void set_high(bool&& _high) { high = _high; }

    // void set_head_node(N* _head) {
    //     job = _head->get_job();
    //     weight = _head->get_weight();
    //     head_node = _head;
    // }

    void set_node_id(NodeId _id) { node_id = _id; }

    void reset() {
        f = DBL_MAX;
        prev = nullptr;
        high = false;
    }

    T get_f() const { return f; }

    T& get_f() { return f; }

    Label<N, T>* get_previous() { return prev; }

    bool get_high() { return high; }

    Job* get_job() { return job; }

    NodeId& get_node_id() { return node_id; }

    Job* get_previous_job() {
        return get_previous() == nullptr ? nullptr : get_previous()->get_job();
    }

    void update_solution(T _f, Label<N, T>* _prev, bool _high = false) {
        f = _f;
        prev = _prev;
        high = _high;
    }

    void update_solution(T _f, Label<N, T>& _node) {
        f = _f;
        prev = _node.prev;
        high = _node.high;
    }

    void update_solution(Label<N, T>& _node) {
        f = _node.f;
        prev = _node.prev;
        high = _node.high;
    }

    N* get_node() const { return head_node; }

    int get_weight() { return weight; }

    void update_label(Label<N, T>* _n, T _f = 0, bool _high = false) {
        if (_high) {
            f = _f;
            prev_job = get_job();
        } else {
            // f = _n->f;
            f = _f;
            prev_job = _n->prev_job;
        }

        high = _high;
        prev = _n;
    }

    Job* get_prev_job() { return prev_job; }
};

#endif  // LABEL_HPP
