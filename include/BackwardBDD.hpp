#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
#include <fmt/core.h>
#include <limits>
#include <range/v3/action/remove_if.hpp>
#include <span>
#include "NodeBdd.hpp"
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"

template <typename T = double>
class BackwardBddBase : public Eval<NodeBdd<T>, OptimalSolution<T>> {
    double* pi{};

   public:
    BackwardBddBase() = default;

    void set_pi(double* _pi) { pi = _pi; }

    [[nodiscard]] const double* get_pi() const { return pi; }

    OptimalSolution<T> get_objective(NodeBdd<T>& n) const {
        OptimalSolution<T> sol(0.0);
        auto*              aux_label = &(n.backward_label[0]);
        auto* table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();

        do {
            auto tmp_node_id = aux_label->get_node_id();
            if (aux_label->get_high()) {
                auto& node = table_tmp->node(tmp_node_id);
                sol.push_job_back(node.get_job(), node.reduced_cost[1]);
            }

            aux_label = aux_label->get_previous();
        } while (aux_label);

        return sol;
    }

    virtual void initializenode(NodeBdd<T>& n) const = 0;
    virtual void initializerootnode(NodeBdd<T>& n) const = 0;
    virtual void evalNode(NodeBdd<T>& n) const = 0;

    BackwardBddBase<T>(const BackwardBddBase<T>&) = default;
    BackwardBddBase<T>(BackwardBddBase<T>&&) noexcept = default;
    BackwardBddBase<T>& operator=(const BackwardBddBase<T>&) = default;
    BackwardBddBase<T>& operator=(BackwardBddBase<T>&&) noexcept = default;
    virtual ~BackwardBddBase<T>() = default;
};

template <typename T = double>
class BackwardBddSimple : public BackwardBddBase<T> {
   public:
    BackwardBddSimple() = default;

    void evalNode(NodeBdd<T>& n) const override {
        auto  table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();
        auto& p0_tmp = table_tmp->node(n[0]);
        auto& p1_tmp = table_tmp->node(n[1]);

        n.reset_reduced_costs();

        const double* dual = BackwardBddBase<T>::get_pi();

        for (auto& list : n.coeff_list) {
            list |= ranges::actions::remove_if([&](auto it) {
                auto aux = it.lock();
                if (aux) {
                    n.adjust_reduced_costs(
                        aux->get_coeff() * dual[aux->get_row()],
                        aux->get_high());
                    return false;
                } else {
                    return true;
                }
            });
        }

        auto obj0 = p0_tmp.backward_label[0].get_f() + n.reduced_cost[0];
        auto obj1 = p1_tmp.backward_label[0].get_f() + n.reduced_cost[1];

        if (obj0 > obj1) {
            n.backward_label[0].backward_update(&(p1_tmp.backward_label[0]),
                                                obj1, true);
        } else {
            n.backward_label[0].backward_update(&(p0_tmp.backward_label[0]),
                                                obj0, false);
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].reset();
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.backward_label[0].get_f() = 0;
    }

    BackwardBddSimple<T>(const BackwardBddSimple<T>&) = default;
    BackwardBddSimple<T>(BackwardBddSimple<T>&&) noexcept = default;
    BackwardBddSimple<T>& operator=(const BackwardBddSimple<T>&) = default;
    BackwardBddSimple<T>& operator=(BackwardBddSimple<T>&&) noexcept = default;
    virtual ~BackwardBddSimple<T>() = default;
};

template <typename T = double>
class BackwardBddCycle : public BackwardBddBase<T> {
   public:
    BackwardBddCycle<T>() = default;

    void evalNode(NodeBdd<T>& n) const override {
        auto* tmp_j = n.get_job();
        auto* table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();
        auto& p0_tmp = table_tmp->node(n[0]);
        auto& p1_tmp = table_tmp->node(n[1]);
        const auto dual = BackwardBddBase<T>::get_pi();

        n.reset_reduced_costs();

        n.adjust_reduced_costs(dual[tmp_j->job], true);
        n.adjust_reduced_costs(0.0, false);

        // for (auto& list : n.coeff_list) {
        //     list |= ranges::actions::remove_if([&](auto it) {
        //         auto aux = it.lock();
        //         if (aux) {
        //             n.adjust_reduced_costs(
        //                 aux->get_coeff() * dual[aux->get_row()],
        //                 aux->get_high());
        //             return false;
        //         } else {
        //             return true;
        //         }
        //     });
        // }

        auto prev_job{p1_tmp.backward_label[0].prev_job_backward()};

        n.backward_label[0].backward_update(
            &(p0_tmp.backward_label[0]),
            p0_tmp.backward_label[0].get_f() + n.reduced_cost[0], false);
        n.backward_label[1].backward_update(
            &(p0_tmp.backward_label[1]),
            p0_tmp.backward_label[1].get_f() + n.reduced_cost[0], false);

        if (prev_job != tmp_j) {
            auto obj1{p1_tmp.backward_label[0].get_f() + n.reduced_cost[1]};
            auto obj2{p1_tmp.backward_label[1].get_f() + n.reduced_cost[1]};

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].prev_job_backward()) {
                    n.backward_label[1].backward_update(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() + n.reduced_cost[0],
                        false);
                }

                n.backward_label[0].backward_update(&(p1_tmp.backward_label[0]),
                                                    obj1, true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].prev_job_backward()) {
                n.backward_label[1].backward_update(&(p1_tmp.backward_label[0]),
                                                    obj1, true);
            } else if (obj2 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].prev_job_backward()) {
                n.backward_label[1].backward_update(&(p1_tmp.backward_label[1]),
                                                    obj2, true);
            }
        } else {
            auto obj1 = p1_tmp.backward_label[1].get_f() + n.reduced_cost[1];

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].prev_job_backward()) {
                    n.backward_label[1].backward_update(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() + n.reduced_cost[0],
                        false);
                }

                n.backward_label[0].backward_update(&(p1_tmp.backward_label[1]),
                                                    obj1, true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].prev_job_backward()) {
                n.backward_label[1].backward_update(&(p1_tmp.backward_label[1]),
                                                    obj1, true);
            }
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].reset();
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.backward_label[0].get_f() = 0.0;
    }

    BackwardBddCycle<T>(const BackwardBddCycle<T>&) = default;
    BackwardBddCycle<T>(BackwardBddCycle<T>&&) noexcept = default;
    BackwardBddCycle<T>& operator=(const BackwardBddCycle<T>&) = default;
    BackwardBddCycle<T>& operator=(BackwardBddCycle<T>&&) noexcept = default;
    virtual ~BackwardBddCycle<T>() = default;
};

#endif  // BACKWARD_BDD_HPP
