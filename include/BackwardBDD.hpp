#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
#include <fmt/core.h>
#include <limits>
#include <span>
#include "NodeBdd.hpp"
#include "NodeBddEval.hpp"
#include "OptimalSolution.hpp"

template <typename T = double>
class BackwardBddBase : public Eval<NodeBdd<T>, OptimalSolution<T>> {
    double* pi{};

   public:
    BackwardBddBase() = default;

    // BackwardBddBase(const BackwardBddBase<T>& src) {}

    void set_pi(double* _pi) { pi = _pi; }

    [[nodiscard]] const double* get_pi() const { return pi; }

    OptimalSolution<T> get_objective(NodeBdd<T>& n) const {
        OptimalSolution<T> sol(0.0);
        auto*              aux_label = &(n.backward_label[0]);
        // auto                  tmp_node_id = aux_label->get_node_id();
        auto table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();

        do {
            auto tmp_node_id = aux_label->get_node_id();
            if (aux_label->get_high()) {
                auto& node = table_tmp->node(tmp_node_id);
                Job*  aux_job = node.get_job();
                sol.push_job_back(aux_job, node.reduced_cost[1]);
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
        auto& p0_tmp = table_tmp->node(n.branch[0]);
        auto& p1_tmp = table_tmp->node(n.branch[1]);

        n.reset_reduced_costs();

        const double* dual = BackwardBddBase<T>::get_pi();

        auto func = [&](auto it) {
            auto aux = it.lock();
            if (aux) {
                n.adjust_reduced_costs(aux->get_coeff() * dual[aux->get_row()],
                                       aux->get_high());
                return false;
            } else {
                return true;
            }
        };

        for (int k = 0; k < 2; k++) {
            auto it = std::remove_if(n.coeff_list[k].begin(),
                                     n.coeff_list[k].end(), func);
            // for (auto it = n.coeff_list[k].begin(); it !=
            // n.coeff_list[k].end();
            //      it++) {
            //     auto aux = it->lock();
            //     if (aux) {
            //         n.adjust_reduced_costs(
            //             aux->get_coeff() * dual[aux->get_row()],
            //             aux->get_high());
            //     }
            // }
            n.coeff_list[k].erase(it, n.coeff_list[k].end());
        }

        T obj0 = p0_tmp.backward_label[0].get_f() + n.reduced_cost[0];
        T obj1 = p1_tmp.backward_label[0].get_f() + n.reduced_cost[1];

        if (obj0 > obj1) {
            n.backward_label[0].update_label(&(p1_tmp.backward_label[0]), obj1,
                                             true);
        } else {
            n.backward_label[0].update_label(&(p0_tmp.backward_label[0]), obj0,
                                             false);
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].update_solution(
            std::numeric_limits<double>::max() / 2, nullptr, false);
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
        auto  tmp_j = n.get_job();
        auto  table_tmp = Eval<NodeBdd<T>, OptimalSolution<T>>::get_table();
        auto& p0_tmp = table_tmp->node(n.branch[0]);
        auto& p1_tmp = table_tmp->node(n.branch[1]);
        const auto dual = BackwardBddBase<T>::get_pi();

        n.reset_reduced_costs();
        auto func = [&](auto it) {
            auto aux = it.lock();
            if (aux) {
                n.adjust_reduced_costs(aux->get_coeff() * dual[aux->get_row()],
                                       aux->get_high());
                return false;
            } else {
                return true;
            }
        };

        for (auto k = 0; k < 2; k++) {
            auto it = std::remove_if(n.coeff_list[k].begin(),
                                     n.coeff_list[k].end(), func);
            n.coeff_list[k].erase(it, n.coeff_list[k].end());
        }

        auto prev_job{p1_tmp.backward_label[0].get_prev_job()};

        n.backward_label[0].update_label(
            &(p0_tmp.backward_label[0]),
            p0_tmp.backward_label[0].get_f() + n.reduced_cost[0], false);
        n.backward_label[1].update_label(
            &(p0_tmp.backward_label[1]),
            p0_tmp.backward_label[1].get_f() + n.reduced_cost[0], false);
        // bool diff = bool_diff_Fij(weight, prev_job, tmp_j);
        // bool diff1 =
        //     bool_diff_Fij(weight, p1_tmp.backward_label[0].get_prev_job(),
        //     tmp_j);

        if (prev_job != tmp_j) {
            T obj1{p1_tmp.backward_label[0].get_f() + n.reduced_cost[1]};
            T obj2{p1_tmp.backward_label[1].get_f() + n.reduced_cost[1]};

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].get_prev_job()) {
                    n.backward_label[1].update_label(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() + n.reduced_cost[0],
                        false);
                }

                n.backward_label[0].update_label(&(p1_tmp.backward_label[0]),
                                                 obj1, true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].get_prev_job()) {
                n.backward_label[1].update_label(&(p1_tmp.backward_label[0]),
                                                 obj1, true);
            } else if (obj2 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].get_prev_job()) {
                n.backward_label[1].update_label(&(p1_tmp.backward_label[1]),
                                                 obj2, true);
            }
        } else {
            T obj1 = p1_tmp.backward_label[1].get_f() + n.reduced_cost[1];

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].get_prev_job()) {
                    n.backward_label[1].update_label(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() + n.reduced_cost[0],
                        false);
                }

                n.backward_label[0].update_label(&(p1_tmp.backward_label[1]),
                                                 obj1, true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].get_prev_job()) {
                n.backward_label[1].update_label(&(p1_tmp.backward_label[1]),
                                                 obj1, true);
            }
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].update_solution(
            std::numeric_limits<double>::max() / 2, nullptr, false);
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
