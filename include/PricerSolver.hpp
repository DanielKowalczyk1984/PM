#ifndef INCLUDE_PRICERSOLVER_HPP
#define INCLUDE_PRICERSOLVER_HPP

#include <iostream>
#include <vector>
#include <PricerEvaluate.hpp>
#include <PricerConstruct.hpp>
#include <tdzdd/DdStructure.hpp>
#include <tdzdd/op/Lookahead.hpp>

struct PricerSolver {
   public:
    tdzdd::DdStructure<2> *                   zdd;
    GPtrArray *                               interval_list;
    int                                       njobs;
    int **sum_p;
    int nlayers;
    tdzdd::DataTable<PricerWeightZDD<double>> zdd_table;
    tdzdd::DataTable<PricerWeightBDD<double>> dd_table;
    tdzdd::DataTable<PricerFarkasZDD<double>> farkas_table;

    PricerSolver(GPtrArray *_interval_list, int &njobs)
        : interval_list(_interval_list),
          njobs(njobs) {
            nlayers = (int) interval_list->len;

            PricerConstruct ps(interval_list);
            zdd = new tdzdd::DdStructure<2>(ps);
            long tmp_size = zdd->size();
            zdd->zddReduce();
            printf("size BDD = %lu, size ZDD= %lu\n",tmp_size, zdd->size());
            init_zdd_table();
            // init_table_farkas();
    };

    PricerSolver(const PricerSolver &other) {
        zdd = new tdzdd::DdStructure<2>;
        *zdd = *(other.zdd);
        interval_list = other.interval_list;
        njobs = other.njobs;
        nlayers = other.nlayers;
        zdd_table.init();
        dd_table.init();
        farkas_table.init();
    }

    void init_tables() {
        init_zdd_table();
        init_table_farkas();
    }

    ~PricerSolver() {
        delete zdd;
    }

    void create_dot_zdd(const char *name) {
        std::ofstream file;
        file.open(name);
        zdd->dumpDot(file);
        file.close();
    }

    void init_table_farkas() {
        tdzdd::NodeTableHandler<2> &node_handler_farkas = zdd->getDiagram();
        farkas_table.init(njobs + 1);
        size_t const m = node_handler_farkas.privateEntity()[0].size();
        farkas_table[0].resize(m);

        for (size_t i = 0; i < m; ++i) {
            farkas_table[0][i].init_terminal_node(i);
        }

        for (int i = 1; i <= njobs; ++i) {
            tdzdd::MyVector<tdzdd::Node<2>> const &node =
                node_handler_farkas.privateEntity()[i];
            size_t const mm = node.size();
            farkas_table[i].resize(mm);

            for (size_t j = 0; j < mm; ++j) {
                farkas_table[i][j].init_node();
            }
        }
    }

    void init_zdd_table() {
        tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
        tdzdd::NodeId &             root = zdd->root();
        zdd_table.init(root.row() + 1);

        /** init table */
        for (int i = root.row(); i >= 0; i--) {
            tdzdd::MyVector<tdzdd::Node<2>> const &node =
                handler.privateEntity()[i];
            size_t const m = node.size();
            zdd_table[i].resize(m);
        }

        /** init root */
        zdd_table[root.row()][root.col()].add_weight(0, 0);

        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();
            int          layer= nlayers - i;
            job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(interval_list,layer);
            Job *job = tmp_pair->j;

            for (size_t j = 0; j < m; j++) {
                tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
                tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

                // for (auto &it : zdd_table[i][j].info_node) {
                //     zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight(it.first);
                //     zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight(it.first
                //     + p[cur_job]);
                // }

                for (my_iterator<double> it = zdd_table[i][j].list.begin();
                     it != zdd_table[i][j].list.end(); it++) {
                    (*it)->n = zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight((*it)->weight,nlayers - cur_node_0.row());
                    (*it)->y = zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight((*it)->weight + job->processingime,nlayers - cur_node_1.row());
                }
            }
        }

        /** init terminal nodes */
        size_t const mm = handler.privateEntity()[0].size();

        for (size_t j = 0; j < mm; j++) {
            zdd_table[0][j].init_terminal_node(j);
        }
    }

    void init_zdd_one_conflict(int v1, int v2, int same) {
        int  ecount_same = 0;
        int  ecount_diff = 0;
        int *elist_same = (int *)NULL;
        int *elist_differ = (int *)NULL;

        if (same) {
            ecount_same = 1;
            elist_same = new int[2];
            elist_same[0] = v1;
            elist_same[1] = v2;
        } else {
            ecount_diff = 1;
            elist_differ = new int[2];
            elist_differ[0] = v1;
            elist_differ[1] = v2;
        }

        ConflictConstraints conflict(nlayers, elist_same, ecount_same,
                                     elist_differ, ecount_diff);
        zdd->zddSubset(tdzdd::ZddLookahead<ConflictConstraints>(conflict));
        zdd->zddReduce();
        delete[] elist_same;
        delete[] elist_differ;
    }


    void init_zdd_conflict_solver(int *elist_same,
                                  int  ecount_same,
                                  int *elist_differ,
                                  int  ecount_differ) {
        if (ecount_same + ecount_differ > 0) {
            // zdd = new tdzdd::DdStructure<2>;
            ConflictConstraints conflict(njobs, elist_same, ecount_same,
                                         elist_differ, ecount_differ);
            zdd->zddSubset(conflict);
            zdd->zddReduce();
        } else {
        }

        init_zdd_table();
        // init_table_farkas();
    }

    void free_zdd_solver(int ecount_same, int ecount_differ) {
        if (ecount_same + ecount_differ > 0) {
            delete zdd;
        }

        zdd_table.init();
        farkas_table.init();
    }

    class Optimal_Solution<double>
    solve_duration_bdd_double(double *pi) {
        return Optimal_Solution<double>();
    }

    class Optimal_Solution<double>
    solve_duration_zdd_double(double *pi) {
        return Optimal_Solution<double>();
    }

    class Optimal_Solution<double>
    solve_weight_bdd_double(double *pi) {
        return Optimal_Solution<double>();
    }

    class Optimal_Solution<double>
    solve_weight_zdd_double(double *pi) {
        return zdd->evaluate_weight(WeightZDDdouble(pi, interval_list, njobs), zdd_table);
    }

    class Optimal_Solution<double>
    solve_farkas_double(double *pi) {
        return Optimal_Solution<double>();
    }

    PricerSolver &
    operator=(PricerSolver const &other) {
        if (this != &other) {
            zdd = new tdzdd::DdStructure<2>;
            *zdd = *(other.zdd);
            interval_list = other.interval_list;
            nlayers = other.nlayers;
            njobs = other.njobs;
            zdd_table = tdzdd::DataTable<PricerWeightZDD<double>>();
            dd_table = tdzdd::DataTable<PricerWeightBDD<double>>();
            farkas_table = tdzdd::DataTable<PricerFarkasZDD<double>>();
        }

        return *this;
    }

    void set_release_due_time(Job *_jobarray) { ; }

    void iterate_zdd() {
        tdzdd::DdStructure<2>::const_iterator it = zdd->begin();

        for (; it != zdd->end(); ++it) {
            std::set<int>::const_iterator i = (*it).begin();

            for (; i != (*it).end(); ++i) {
                std::cout << njobs - *i << " ";
            }

            std::cout << std::endl;
        }
    }
};

#endif  // INCLUDE_PRICERSOLVER_HPP
