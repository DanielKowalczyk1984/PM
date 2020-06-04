#ifndef LABEL_HPP
#define LABEL_HPP

#include "OptimalSolution.hpp"

template<typename N, typename T>
class Label
{
    private:
        Label<N,T>* prev;
        bool high;
        N* head_node;

    public:
        T f;
        Job* prev_job;
        /**
         * Constructor
         */
        Label(T& _f, Label<N,T>*& _prev, bool& _high) :
            prev(_prev),
            high(_high),
            head_node(nullptr),
            f(_f),
            prev_job(nullptr) {};

        Label() :
            prev(nullptr),
            high(false),
            head_node(nullptr),
            f(DBL_MAX),
            prev_job(nullptr) {};

        explicit Label(N* _head_node) :
            prev(nullptr),
            high(false),
            head_node(_head_node),
            f(DBL_MAX),
            prev_job(nullptr) {};

        Label<N,T>(const Label<N,T>& src) = default;
        Label<N,T>(Label<N,T>&& src) = default;
        Label<N,T>& operator=(const Label<N,T>& src) = default;
        Label<N,T>& operator=(Label<N,T>&& src) = default;

        void set_previous(Label<N,T>*&& _prev)
        {
            prev = _prev;
        }

        void set_f(T _f)
        {
            f = _f;
        }

        void set_high(bool&& _high)
        {
            high = _high;
        }

        void set_head_node(N* _head)
        {
            head_node = _head;
        }

        void reset()
        {
            f = DBL_MAX;
            prev = nullptr;
            high = false;
        }

        T get_f() const
        {
            return f;
        }

        Label<N,T>* get_previous()
        {
            return prev;
        }

        bool get_high()
        {
            return high;
        }

        Job* get_job()
        {
            return head_node->get_job();
        }

        Job* get_previous_job()
        {
            return get_previous() == nullptr ? nullptr : get_previous()->get_job();
        }

        void update_solution(T _f, Label<N,T>*&& _prev, bool&& _high)
        {
            f = _f;
            prev = _prev;
            high = _high;
        }

        void update_solution(Label<N,T>& _node)
        {
            f = _node.f;
            prev = _node.prev;
            high = _node.high;
        }

        N* get_node() const
        {
            return head_node;
        }

        int get_weight()
        {
            return head_node->get_weight();
        }

        void update_label(Label<N,T>* _n, T _f = 0, bool _high = false)
        {
            if (_high) {
                f = _f;
                prev_job = get_job();
            } else {
                f = _n->f;
                prev_job = _n->prev_job;
            }

            high = _high;
            prev = _n;
        }

        Job* get_prev_job()
        {
            return prev_job;
        }
};


#endif // LABEL_HPP


