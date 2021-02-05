#ifndef FORWARD_BDD_HPP
#define FORWARD_BDD_HPP
// #include <tdzdd/DdEval.hpp>
#include <fmt/core.h>
#include <NodeBdd.hpp>
#include <OptimalSolution.hpp>
#include <array>
#include "NodeBddEval.hpp"
#include "math.h"

template <typename T = double>
class ForwardBddBase : public Eval<NodeBdd<T>, OptimalSolution<T>> {
    double* pi{};

   public:
    ForwardBddBase() = default;

    ForwardBddBase<T>(const ForwardBddBase<T>& src) = default;
    ForwardBddBase<T>(ForwardBddBase<T>&&) noexcept = default;
    ForwardBddBase<T>& operator=(ForwardBddBase<T>&&) noexcept = default;
    ForwardBddBase<T>& operator=(const ForwardBddBase<T>&) = default;
    ~ForwardBddBase<T>() = default;

    void set_pi(double* _pi) { pi = _pi; }

    [[nodiscard]] const double* get_pi() const { return pi; }

    virtual void initializenode(NodeBdd<T>& n) const = 0;

    virtual void initializerootnode(NodeBdd<T>& n) const = 0;

    virtual void evalNode(NodeBdd<T>& n) const = 0;

    OptimalSolution<T> get_objective(NodeBdd<T>& n) const {
        OptimalSolution<T>    sol(0);
        Label<NodeBdd<T>, T>* ptr_node = &(n.forward_label[0]);
        auto table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();

        while (ptr_node->get_previous() != nullptr) {
            Label<NodeBdd<T>, T>* aux_prev_node = ptr_node->get_previous();
            auto& node = table_tmp->node(aux_prev_node->get_node_id());
            Job*  aux_job = node.get_job();
            sol.C_max += aux_job->processing_time;
            sol.push_job_back(aux_job, node.get_weight(), node.reduced_cost[1]);
            ptr_node = aux_prev_node;
        }

        sol.reverse_jobs();

        return sol;
    }
};

template <typename T = double>
class ForwardBddCycle : public ForwardBddBase<T> {
   public:
    ForwardBddCycle() = default;

    ForwardBddCycle<T>(const ForwardBddCycle<T>& src) = default;
    ForwardBddCycle<T>(ForwardBddCycle<T>&&) noexcept = default;
    ForwardBddCycle<T>& operator=(const ForwardBddCycle<T>&) = default;
    ForwardBddCycle<T>& operator=(ForwardBddCycle<T>&&) noexcept = default;
    ~ForwardBddCycle<T>() = default;

    void initializenode(NodeBdd<T>& n) const override {
        if (n.get_weight() == 0) {
            n.forward_label[0].update_solution(0, nullptr, false);
            n.forward_label[1].update_solution(DBL_MAX / 2, nullptr, false);
        } else {
            n.forward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
            n.forward_label[1].update_solution(DBL_MAX / 2, nullptr, false);
        }
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.forward_label[0].get_f() = 0;
        n.forward_label[1].set_f(DBL_MAX / 2);
    }

    void evalNode(NodeBdd<T>& n) const override {
        Job* tmp_j = n.get_job();
        assert(tmp_j != nullptr);
        double result = NAN;
        bool   diff = false;

        int   weight = n.get_weight();
        T     g;
        auto  table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();
        auto& p0 = table_tmp->node(n.branch[0]);
        auto& p1 = table_tmp->node(n.branch[1]);
        n.reset_reduced_costs();
        const double* dual = ForwardBddBase<T>::get_pi();

        for (int k = 0; k < 2; k++) {
            for (auto it = n.coeff_list[k].begin(); it != n.coeff_list[k].end();
                 it++) {
                auto aux = it->lock();
                if (aux) {
                    n.adjust_reduced_costs(
                        aux->get_coeff() * dual[aux->get_row()],
                        aux->get_high());
                }
            }
        }
        result = n.reduced_cost[1];

        /**
         * High edge calculation
         */
        Job* prev = n.forward_label[0].get_previous_job();
        Job* aux1 = p1.forward_label[0].get_previous_job();
        // // diff = (prev == nullptr) ? true
        //                          : (value_diff_Fij(weight, tmp_j, prev) >=
        //                          0);

        if (prev != tmp_j
            // && diff
        ) {
            g = n.forward_label[0].get_f() + result;
            if (g < p1.forward_label[0].get_f()) {
                if (aux1 != tmp_j) {
                    p1.forward_label[1].update_solution(p1.forward_label[0]);
                }
                p1.forward_label[0].update_solution(g, &(n.forward_label[0]),
                                                    true);
            } else if ((g < p1.forward_label[1].get_f()) && (aux1 != tmp_j)) {
                p1.forward_label[1].update_solution(g, &(n.forward_label[0]),
                                                    true);
            }
        } else {
            g = n.forward_label[1].get_f() + result;
            prev = n.forward_label[1].get_previous_job();
            // diff = (prev == nullptr)
            //    ? true
            //    : (value_diff_Fij(weight, tmp_j, prev) >= 0);

            // if (diff) {
            if (g < p1.forward_label[0].get_f()) {
                if (aux1 != tmp_j) {
                    p1.forward_label[1].update_solution(p1.forward_label[0]);
                }
                p1.forward_label[0].update_solution(g, &(n.forward_label[1]),
                                                    true);
            } else if ((g < p1.forward_label[1].get_f()) && (aux1 != tmp_j)) {
                p1.forward_label[1].update_solution(g, &(n.forward_label[1]),
                                                    true);
            }
            // }
        }

        /**
         * Low edge calculation
         */
        aux1 = p0.forward_label[0].get_previous_job();
        g = n.forward_label[0].get_f() + n.reduced_cost[0];
        auto g1 = n.forward_label[1].get_f() + n.reduced_cost[0];
        if (g < p0.forward_label[0].get_f()) {
            if (prev != aux1) {
                p0.forward_label[1].update_solution(p0.forward_label[0]);
            }
            p0.forward_label[0].update_solution(g, n.forward_label[0]);
            if (g1 < p0.forward_label[1].get_f()) {
                p0.forward_label[1].update_solution(g1, n.forward_label[1]);
            }
        } else if ((g < p0.forward_label[1].get_f()) && (aux1 != prev)) {
            p0.forward_label[1].update_solution(g, n.forward_label[0]);
        } else if ((g1 < p0.forward_label[1].get_f())) {
            p0.forward_label[1].update_solution(g1, n.forward_label[1]);
        }
    }
};

template <typename T = double>
class ForwardBddSimple : public ForwardBddBase<T> {
   public:
    ForwardBddSimple<T>() = default;

    ForwardBddSimple(const ForwardBddSimple<T>&) = default;
    ForwardBddSimple<T>& operator=(const ForwardBddSimple<T>&) = default;
    ForwardBddSimple(ForwardBddSimple<T>&&) noexcept = default;
    ForwardBddSimple<T>& operator=(ForwardBddSimple<T>&&) noexcept = default;
    ~ForwardBddSimple<T>() = default;

    void initializenode(NodeBdd<T>& n) const override {
        if (n.get_weight() == 0) {
            n.forward_label[0].update_solution(0, nullptr, false);
        } else {
            n.forward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
        }
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.forward_label[0].get_f() = 0;
    }

    void evalNode(NodeBdd<T>& n) const override {
        T     g;
        auto  table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();
        auto& p0 = table_tmp->node(n.branch[0]);
        auto& p1 = table_tmp->node(n.branch[1]);
        n.reset_reduced_costs();
        const double* dual = ForwardBddBase<T>::get_pi();

        for (int k = 0; k < 2; k++) {
            for (auto it = n.coeff_list[k].begin(); it != n.coeff_list[k].end();
                 it++) {
                auto aux = it->lock();
                if (aux) {
                    n.adjust_reduced_costs(
                        aux->get_coeff() * dual[aux->get_row()],
                        aux->get_high());
                }
            }
        }

        /**
         * High edge calculation
         */
        g = n.forward_label[0].get_f() + n.reduced_cost[1];
        if (g < p1.forward_label[0].get_f()) {
            p1.forward_label[0].update_solution(g, &(n.forward_label[0]), true);
        }

        /**
         * Low edge calculation
         */
        g = n.forward_label[0].get_f() + n.reduced_cost[0];
        double result = g;
        if (g < p0.forward_label[0].get_f()) {
            p0.forward_label[0].update_solution(result, n.forward_label[0]);
        }
    }
};

#endif  // FORWARD_BDD_HPP
