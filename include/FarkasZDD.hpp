#include <solution.h>
#include <interval.h>
#include <tdzdd/DdEval.hpp>
#include <cfloat>
#include <vector>
#include <OptimalSolution.hpp>
#include <tdzdd/dd/NodeTable.hpp>

template <typename T>
class PricerFarkasZDD {
public:
    T    obj;
    bool take;

    PricerFarkasZDD() : obj(0), take(0) {};

    ~PricerFarkasZDD() {};

    void init_terminal_node(int one) { obj = one ? 0.0 : 1871286761.0; }

    void init_node() { take = false; }
};

/**
 * Farkas
 */
template <typename E, typename T>
class FarkasZDD
    : public tdzdd::DdEval<E, PricerFarkasZDD<T>, Optimal_Solution<T>> {
    T   *pi;
    GPtrArray *interval_list;
    int nlayers;
    int nbjobs;
    job_interval_pair *tmp_pair;
    Job *tmp_j;
    interval *tmp_interval;

public:
    FarkasZDD(T *_pi, GPtrArray *_interval_list, int _nbjobs)
        : pi(_pi),
          interval_list(_interval_list),
          nbjobs(_nbjobs) {};

    void evalTerminal(PricerFarkasZDD<T>& n) { n.obj = pi[nbjobs]; }

    void evalNode(PricerFarkasZDD<T> *n,
                  int                 i,
                  tdzdd::DdValues<PricerFarkasZDD<T>, 2>& values) const
    {
        int j = nlayers - i;
        assert(j >= 0 && j <= nbjobs - 1);
        tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list, j);
        tmp_j = tmp_pair->j;
        PricerFarkasZDD<T> *n0 = values.get_ptr(0);
        PricerFarkasZDD<T> *n1 = values.get_ptr(1);

        if (n0->obj < n1->obj + pi[tmp_j->job]) {
            n->obj = n0->obj;
            n->take = false;
        } else {
            n->obj = n1->obj + pi[tmp_j->job];
            n->take = true;
        }
    }

    void initializenode(PricerFarkasZDD<T>& n) { n.take = false; }

    Optimal_Solution<T> get_objective(
        tdzdd::NodeTableHandler<2>            diagram,
        tdzdd::DataTable<PricerFarkasZDD<T>> *data_table,
        const tdzdd::NodeId                  *f)
    {
        Optimal_Solution<T> sol;
        sol.obj = (*data_table)[f->row()][f->col()].obj;
        sol.jobs = g_ptr_array_new();
        tdzdd::NodeId cur_node = *f;
        int           j = nlayers - cur_node.row();

        while (cur_node.row() != 0) {
            if ((*data_table)[cur_node.row()][cur_node.col()].take) {
                g_ptr_array_add(sol.jobs, tmp_j);
                sol.C_max += tmp_j->processingime;
                sol.cost += value_Fj(sol.C_max, tmp_j);
                cur_node = diagram.privateEntity().child(cur_node, 1);
            } else {
                cur_node = diagram.privateEntity().child(cur_node, 0);
            }

            j = nlayers - cur_node.row();
        }

        return sol;
    }
};
