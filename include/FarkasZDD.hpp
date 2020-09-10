#include <BackwardBDD.hpp>

template <typename T = double>
class BackwardBddFarkas : public BackwardBddBase<T> {
   public:
    BackwardBddFarkas() : BackwardBddBase<T>(){};

    void evalNode(NodeBdd<T>& n) const override {
        n.reset_reduced_costs_farkas();

        const double* dual = BackwardBddBase<T>::get_pi();
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

        // Job *tmp_j = n.get_job();
        NodeBdd<T>* p0 = n.child[0];
        NodeBdd<T>* p1 = n.child[1];
        // T result = -pi[tmp_j->job];

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
        n.backward_label[0].f = 0.0;
    }

    OptimalSolution<T> get_objective(NodeBdd<T>& n) const {
        OptimalSolution<T> sol(0.0);

        NodeBdd<T>* aux_node = &n;
        Job*        aux_job = n.get_job();

        while (aux_job) {
            if (aux_node->backward_label[0].get_high()) {
                sol.push_job_back_farkas(aux_job, n.reduced_cost[1]);
                aux_node = aux_node->child[1];
                aux_job = aux_node->get_job();
            } else {
                aux_node = aux_node->child[0];
                aux_job = aux_node->get_job();
            }
        }

        return sol;
    }
};
/**
 * Farkas
 */