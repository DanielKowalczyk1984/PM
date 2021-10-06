#ifndef __FARKASZDD_H__
#define __FARKASZDD_H__

#include "BackwardBDD.hpp"           // for BackwardBddBase
#include "ModernDD/NodeBddEval.hpp"  // for Eval
#include "NodeBdd.hpp"               // for NodeBdd
#include "PricingSolution.hpp"       // for PricingSolution

class BackwardBddFarkas : public BackwardBddBase {
   public:
    BackwardBddFarkas() = default;

    void evalNode(NodeBdd& n) const override {
        n.reset_reduced_costs_farkas();

        const double* dual = BackwardBddBase::get_pi();
        for (auto& list : n.get_coeff_list()) {
            for (auto& it : list) {
                auto aux = it.lock();
                if (aux) {
                    n.adjust_reduced_costs(
                        aux->get_coeff() * dual[aux->get_row()],
                        aux->get_high());
                }
            }
        }

        auto  table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto& p0 = table_tmp->node(n[0]);
        auto& p1 = table_tmp->node(n[1]);

        auto obj0 = p0.backward_label[0].get_f() + n.get_reduced_cost()[0];
        auto obj1 = p1.backward_label[0].get_f() + n.get_reduced_cost()[1];

        if (obj0 > obj1) {
            n.backward_label[0].backward_update(&(p1.backward_label[0]), obj1,
                                                true);
        } else {
            n.backward_label[0].backward_update(&(p0.backward_label[0]), obj0,
                                                false);
        }
    }

    void initialize_node(NodeBdd& n) const override {
        n.backward_label[0].reset();
    }

    void initializerootnode(NodeBdd& n) const override {
        n.backward_label[0].get_f() = 0.0;
    }
};
#endif  // __FARKASZDD_H__