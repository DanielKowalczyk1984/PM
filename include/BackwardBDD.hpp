#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
#include <array>
#include <range/v3/action/remove_if.hpp>  // for remove_if
#include <span>                           //for span
#include "ModernDD/NodeBddEval.hpp"       // for Eval
#include "NodeBdd.hpp"                    // for NodeBdd
#include "PricingSolution.hpp"            // for PricingSolution
#include "BddCoeff.hpp"

class BackwardBddBase : public Eval<NodeBdd, PricingSolution> {
    const double* pi{};

   public:
    BackwardBddBase() = default;

    void set_pi(double* _pi) { pi = _pi; }
    void set_pi(std::span<const double>& _pi) { pi = _pi.data(); }

    [[nodiscard]] const double* get_pi() const { return pi; }

    PricingSolution get_objective(NodeBdd& n) const override {
        PricingSolution sol(0.0);
        auto*           aux_label = &(n.backward_label[0]);
        auto*           table_tmp = Eval<NodeBdd, PricingSolution>::get_table();

        do {
            auto& tmp_node_id = aux_label->get_node_id();
            if (aux_label->get_high()) {
                auto& node = table_tmp->node(tmp_node_id);
                sol.push_job_back(node.get_job(), node.get_reduced_cost()[1]);
            }

            aux_label = aux_label->get_previous();
        } while (aux_label);

        return sol;
    }

    void initialize_node(NodeBdd& n) const override = 0;
    void initialize_root_node(NodeBdd& n) const override = 0;
    void evalNode(NodeBdd& n) const override = 0;

    BackwardBddBase(const BackwardBddBase&) = default;
    BackwardBddBase(BackwardBddBase&&) noexcept = default;
    BackwardBddBase& operator=(const BackwardBddBase&) = default;
    BackwardBddBase& operator=(BackwardBddBase&&) noexcept = default;
    ~BackwardBddBase() override = default;
};

class BackwardBddSimple : public BackwardBddBase {
   public:
    BackwardBddSimple() = default;

    void evalNode(NodeBdd& n) const override {
        auto  table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto& p0_tmp = table_tmp->node(n[0]);
        auto& p1_tmp = table_tmp->node(n[1]);

        n.reset_reduced_costs();

        const double* dual = BackwardBddBase::get_pi();

        for (auto& list : n.get_coeff_list()) {
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

        auto obj0 = p0_tmp.backward_label[0].get_f() + n.get_reduced_cost()[0];
        auto obj1 = p1_tmp.backward_label[0].get_f() + n.get_reduced_cost()[1];

        if (obj0 > obj1) {
            n.backward_label[0].backward_update(&(p1_tmp.backward_label[0]),
                                                obj1, true);
        } else {
            n.backward_label[0].backward_update(&(p0_tmp.backward_label[0]),
                                                obj0, false);
        }
    }

    void initialize_node(NodeBdd& n) const override {
        n.backward_label[0].reset();
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.backward_label[0].get_f() = 0;
    }

    BackwardBddSimple(const BackwardBddSimple&) = default;
    BackwardBddSimple(BackwardBddSimple&&) noexcept = default;
    BackwardBddSimple& operator=(const BackwardBddSimple&) = default;
    BackwardBddSimple& operator=(BackwardBddSimple&&) noexcept = default;
    ~BackwardBddSimple() override = default;
};

class BackwardBddCycle : public BackwardBddBase {
   public:
    BackwardBddCycle() = default;

    void evalNode(NodeBdd& n) const override {
        auto*      tmp_j = n.get_job();
        auto*      table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto&      p0_tmp = table_tmp->node(n[0]);
        auto&      p1_tmp = table_tmp->node(n[1]);
        const auto dual = BackwardBddBase::get_pi();

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
            p0_tmp.backward_label[0].get_f() + n.get_reduced_cost()[0], false);
        n.backward_label[1].backward_update(
            &(p0_tmp.backward_label[1]),
            p0_tmp.backward_label[1].get_f() + n.get_reduced_cost()[0], false);

        if (prev_job != tmp_j) {
            auto obj1{p1_tmp.backward_label[0].get_f() +
                      n.get_reduced_cost()[1]};
            auto obj2{p1_tmp.backward_label[1].get_f() +
                      n.get_reduced_cost()[1]};

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].prev_job_backward()) {
                    n.backward_label[1].backward_update(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() +
                            n.get_reduced_cost()[0],
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
            auto obj1 =
                p1_tmp.backward_label[1].get_f() + n.get_reduced_cost()[1];

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].prev_job_backward()) {
                    n.backward_label[1].backward_update(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() +
                            n.get_reduced_cost()[0],
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

    void initialize_node(NodeBdd& n) const override {
        n.backward_label[0].reset();
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.backward_label[0].get_f() = 0.0;
    }

    BackwardBddCycle(const BackwardBddCycle&) = default;
    BackwardBddCycle(BackwardBddCycle&&) noexcept = default;
    BackwardBddCycle& operator=(const BackwardBddCycle&) = default;
    BackwardBddCycle& operator=(BackwardBddCycle&&) noexcept = default;
    ~BackwardBddCycle() override = default;
};

#endif  // BACKWARD_BDD_HPP
