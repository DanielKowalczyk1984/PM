#include <limits>
#include "BackwardBDD.hpp"
#include "NodeBddEval.hpp"

template <typename T = double>
class BackwardBddFarkas : public BackwardBddBase<T> {
   public:
    BackwardBddFarkas() = default;

    void evalNode(NodeBdd<T>& n) const override {
        n.reset_reduced_costs_farkas();

        const double* dual = BackwardBddBase<T>::get_pi();
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

        auto  table_tmp = Eval<NodeBdd<T>, PricingSolution<T>>::get_table();
        auto& p0 = table_tmp->node(n[0]);
        auto& p1 = table_tmp->node(n[1]);

        T obj0 = p0.backward_label[0].get_f() + n.get_reduced_cost()[0];
        T obj1 = p1.backward_label[0].get_f() + n.get_reduced_cost()[1];

        if (obj0 > obj1) {
            n.backward_label[0].backward_update(&(p1.backward_label[0]), obj1,
                                                true);
        } else {
            n.backward_label[0].backward_update(&(p0.backward_label[0]), obj0,
                                                false);
        }
    }

    void initializenode(NodeBdd<T>& n) const override {
        n.backward_label[0].reset();
    }

    void initializerootnode(NodeBdd<T>& n) const override {
        n.backward_label[0].get_f() = 0.0;
    }
};
/**
 * Farkas
 */