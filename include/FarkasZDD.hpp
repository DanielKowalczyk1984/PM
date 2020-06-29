#include <BackwardBDD.hpp>

template<typename E, typename T>
class BackwardBddFarkas: public BackwardBddBase<E, T> {
public:
    BackwardBddFarkas() : BackwardBddBase<E, T>() {
    };

    void evalNode(NodeBdd<T> &n) const override {
        // Job *tmp_j = n.get_job();
        NodeBdd<T> *p0 = n.child[0];
        NodeBdd<T> *p1 = n.child[1];
        // T result = -pi[tmp_j->job];

        T obj0 = p0->backward_label[0].get_f() + n.reduced_cost[0];
        T obj1 = p1->backward_label[0].get_f() + n.reduced_cost[1];

        if (obj0 > obj1) {
            n.backward_label[0].update_solution(obj1, nullptr, true);
        } else {
            n.backward_label[0].update_solution(obj0, nullptr, false);
        }

    }

    void initializenode(NodeBdd<T> &n) const override {
        n.backward_label[0].update_solution(DBL_MAX / 2, nullptr, false);
    }

    void initializerootnode(NodeBdd<T> &n) const override {
        n.backward_label[0].f = 0.0;
    }

    OptimalSolution<T> get_objective(NodeBdd<T> &n) const {
        OptimalSolution<T> sol(0.0);

        NodeBdd<T> *aux_node = &n;
        Job *aux_job =  n.get_job();

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

    OptimalSolution<T> getValue(NodeBdd<T> const &n) override {
        OptimalSolution<T> sol;
        return sol;
    }

};
/**
 * Farkas
 */