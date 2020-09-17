#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
#include <fmt/core.h>
#include <NodeBdd.hpp>
#include <OptimalSolution.hpp>
#include "NodeBddEval.hpp"

template <typename T = double>
class BackwardBddBase : public Eval<NodeBdd<T>, OptimalSolution<T>> {
    OriginalModel<>* original_model;
    double*          pi;

   public:
    BackwardBddBase()
        : Eval<NodeBdd<T>, OptimalSolution<T>>(),
          original_model(nullptr),
          pi(nullptr) {}

    // BackwardBddBase(const BackwardBddBase<T>& src) {}

    void set_pi(double* _pi) { pi = _pi; }

    const double* get_pi() const { return pi; }

    virtual void initializenode(NodeBdd<T>& n) const = 0;
    virtual void initializerootnode(NodeBdd<T>& n) const = 0;
    virtual void evalNode(NodeBdd<T>& n) const = 0;
    BackwardBddBase<T>(const BackwardBddBase<T>&) = default;
    BackwardBddBase<T>(BackwardBddBase<T>&&) = default;
    BackwardBddBase<T>& operator=(const BackwardBddBase<T>&) = default;
    BackwardBddBase<T>& operator=(BackwardBddBase<T>&&) = default;
};

template <typename T = double>
class BackwardBddSimple : public BackwardBddBase<T> {
   public:
    BackwardBddSimple() : BackwardBddBase<T>(){};

    void evalNode(NodeBdd<T>& n) const override {
        NodeBdd<T>* p0 = n.child[0];
        NodeBdd<T>* p1 = n.child[1];

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

        T obj0 = p0->backward_label[0].get_f() + n.reduced_cost[0];
        T obj1 = p1->backward_label[0].get_f() + n.reduced_cost[1];

        if (obj0 > obj1) {
            n.backward_label[0].update_solution(obj1, nullptr, true);
        } else {
            n.backward_label[0].update_solution(obj0, nullptr, false);
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.backward_label[0].f = 0;
    }

    OptimalSolution<T> get_objective(NodeBdd<T>& n) const {
        OptimalSolution<T> sol(0.0);

        NodeBdd<T>* aux_node = &n;
        Job*        aux_job = n.get_job();

        while (aux_job) {
            if (aux_node->backward_label[0].get_high()) {
                sol.push_job_back(aux_job, aux_node->reduced_cost[1]);
                aux_node = aux_node->child[1];
                aux_job = aux_node->get_job();
            } else {
                aux_node = aux_node->child[0];
                aux_job = aux_node->get_job();
            }
        }

        return sol;
    }

    BackwardBddSimple<T>(const BackwardBddSimple<T>&) = default;
    BackwardBddSimple<T>(BackwardBddSimple<T>&&) = default;
    BackwardBddSimple<T>& operator=(const BackwardBddSimple<T>&) = default;
    BackwardBddSimple<T>& operator=(BackwardBddSimple<T>&&) = default;
};

template <typename T = double>
class BackwardBddCycle : public BackwardBddBase<T> {
   public:
    BackwardBddCycle() : BackwardBddBase<T>(){};

    void evalNode(NodeBdd<T>& n) const override {
        auto tmp_j = n.get_job();
        // int         weight{n.get_weight()};
        auto       p0{n.child[0]};
        auto       p1{n.child[1]};
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

        auto prev_job{p1->backward_label[0].get_prev_job()};

        n.backward_label[0].update_label(
            &(p0->backward_label[0]),
            p0->backward_label[0].get_f() + n.reduced_cost[0], false);
        n.backward_label[1].update_label(
            &(p0->backward_label[1]),
            p0->backward_label[1].get_f() + n.reduced_cost[0], false);
        // bool diff = bool_diff_Fij(weight, prev_job, tmp_j);
        // bool diff1 =
        //     bool_diff_Fij(weight, p1->backward_label[0].get_prev_job(),
        //     tmp_j);

        if (prev_job != tmp_j) {
            T obj1{p1->backward_label[0].get_f() + n.reduced_cost[1]};
            T obj2{p1->backward_label[1].get_f() + n.reduced_cost[1]};

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].get_prev_job()) {
                    n.backward_label[1].update_label(
                        &(p0->backward_label[0]),
                        p0->backward_label[0].get_f() + n.reduced_cost[0],
                        false);
                }

                n.backward_label[0].update_label(&(p1->backward_label[0]), obj1,
                                                 true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].get_prev_job()) {
                n.backward_label[1].update_label(&(p1->backward_label[0]), obj1,
                                                 true);
            } else if (obj2 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].get_prev_job()) {
                n.backward_label[1].update_label(&(p1->backward_label[1]), obj2,
                                                 true);
            }
        } else {
            T obj1 = p1->backward_label[1].get_f() + n.reduced_cost[1];

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].get_prev_job()) {
                    n.backward_label[1].update_label(
                        &(p0->backward_label[0]),
                        p0->backward_label[0].get_f() + n.reduced_cost[0],
                        false);
                }

                n.backward_label[0].update_label(&(p1->backward_label[1]), obj1,
                                                 true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].get_prev_job()) {
                n.backward_label[1].update_label(&(p1->backward_label[1]), obj1,
                                                 true);
            }
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.backward_label[0].f = 0.0;
    }

    OptimalSolution<T> get_objective(NodeBdd<T>& n) const {
        OptimalSolution<T>    sol(0.0);
        Label<NodeBdd<T>, T>* aux_label = &(n.backward_label[0]);

        while (aux_label) {
            if (aux_label->get_high()) {
                Job* aux_job = aux_label->get_job();
                auto node = aux_label->get_node();
                sol.push_job_back(aux_job, node->reduced_cost[1]);
            }

            aux_label = aux_label->get_previous();
        }

        return sol;
    }
    BackwardBddCycle<T>(const BackwardBddCycle<T>&) = default;
    BackwardBddCycle<T>(BackwardBddCycle<T>&&) = default;
    BackwardBddCycle<T>& operator=(const BackwardBddCycle<T>&) = default;
    BackwardBddCycle<T>& operator=(BackwardBddCycle<T>&&) = default;
};

#endif  // BACKWARD_BDD_HPP
