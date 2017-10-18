#ifndef INCLUDE_PRICERSOLVER_HPP
#define INCLUDE_PRICERSOLVER_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <PricerEvaluate.hpp>
#include <PricerConstruct.hpp>
#include <tdzdd/DdStructure.hpp>
#include <tdzdd/op/Lookahead.hpp>
#include <gurobi_c++.h>
using namespace std;

struct PricerSolver {
public:
    tdzdd::DdStructure<2>                    *zdd;
    GPtrArray                                *interval_list;
    GPtrArray *jobs;
    int                                       njobs;
    int **sum_p;
    int nlayers;
    int remove_edges;
    int remove_nodes;
    tdzdd::DataTable<PricerWeightZDD<double>> zdd_table;
    tdzdd::DataTable<PricerFarkasZDD<double>> farkas_table;
    std::vector<edge<double>> edges;
    GRBEnv *env;
    GRBModel *model;

    PricerSolver(GPtrArray *_interval_list,
                 GPtrArray *_jobs): interval_list(_interval_list), jobs(_jobs)
    {
        njobs = jobs->len;
        nlayers = (int) interval_list->len;
        PricerConstruct ps(interval_list);
        zdd = new tdzdd::DdStructure<2>(ps);
        long tmp_size = zdd->size();
        zdd->zddReduce();
        init_zdd_table();
        printf("size BDD = %lu, size ZDD= %lu\n", tmp_size, zdd->size());
        remove_edges = 0;
        remove_nodes = 0;
        env = new GRBEnv("gurobi.log");
        model = new GRBModel(*env);
    };

    PricerSolver(const PricerSolver& other)
    {
        zdd = new tdzdd::DdStructure<2>;
        *zdd = *(other.zdd);
        interval_list = other.interval_list;
        njobs = other.njobs;
        nlayers = other.nlayers;
        zdd_table.init();
        farkas_table.init();
    }

    void build_mip()
    {
        try {
            /** Initialize MIP model builded with BDD */
            std::vector<edge<double>> e;

            printf("test model\n");
            model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

            /** adding variables */
            for (auto i = edges.begin(); i != edges.end(); ++i) {
                if (i->job) {
                    if(i->out->calc) {
                        i->v = model->addVar(0.0, 1.0, i->cost, GRB_BINARY);
                        e.push_back(*i);
                    } else {
                        i->v = model->addVar(0.0, 0.0, i->cost, GRB_BINARY);
                        e.push_back(*i);
                    }
                } else {
                    if(i->out->calc0) {
                        i->v = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
                        e.push_back(*i);
                    } else {
                        i->v = model->addVar(0.0, 0.0, i->cost, GRB_BINARY);
                        e.push_back(*i);
                    }
                }
            }

            model->update();
            env->set(GRB_IntParam_Threads, 1);
            model->set(GRB_IntParam_Presolve, 2);
            model->set(GRB_IntParam_MIPFocus, 3);
            model->set(GRB_IntParam_Cuts, 2);

            Job *job;

            /** assign the jobs to some path */
            for (unsigned i = 0; i < jobs->len; ++i) {
                job = (Job *) g_ptr_array_index(jobs, i);
                GRBLinExpr assignment = 0;

                for (auto i = e.begin(); i != e.end(); ++i) {
                    if (i->job == job) {
                        assignment += i->v;
                    }
                }

                model->addConstr(assignment, GRB_EQUAL, 1.0);
            }

            tdzdd::NodeId&              root = zdd->root();

            /** Calculate the distance from  the origin to the given node */
            printf("constructing flow constraints\n");
            for (int i = root.row() - 1; i > 0; i--) {
                size_t const m = zdd_table[i].size();

                for (size_t j = 0; j < m; j++) {
                    for (auto it = zdd_table[i][j].list.begin(); it != zdd_table[i][j].list.end(); it++) {
                        GRBLinExpr expr = 0;

                        for (auto v = e.begin(); v != e.end(); ++v) {
                            if (v->out == *it) {
                                expr -= v->v;
                            }

                            if (v->in == *it) {
                                expr += v->v;
                            }
                        }

                        model->addConstr(expr, GRB_EQUAL, 0.0);
                    }
                }
            }

            GRBLinExpr expr = 0;
            for (auto it = zdd_table[0][1].list.begin(); it != zdd_table[0][1].list.end();
                    it++) {
                for (auto v = e.begin(); v != e.end(); v++) {
                    if (v->in == *it) {
                        expr += v->v;
                    }
                }
            }

            model->addConstr(expr, GRB_EQUAL, 4.0);
            printf("begin solving mip\n");
            model->write("test.lp");
            model->optimize();
            // model->computeIIS();
        } catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Exception during optimization" << endl;
        }


    }

    PricerSolver& operator=(PricerSolver const& other)
    {
        if (this != &other) {
            zdd = new tdzdd::DdStructure<2>;
            *zdd = *(other.zdd);
            interval_list = other.interval_list;
            nlayers = other.nlayers;
            njobs = other.njobs;
            zdd_table = tdzdd::DataTable<PricerWeightZDD<double>>();
            farkas_table = tdzdd::DataTable<PricerFarkasZDD<double>>();
        }

        return *this;
    }

    ~PricerSolver() noexcept
    {
        delete zdd;
        if(model) {
            delete model;
        }

        if(env) {
            delete env;
        }
    }

    void init_tables()
    {
        init_zdd_table();
        init_table_farkas();
    }


    void create_dot_zdd(const char *name)
    {
        std::ofstream file;
        file.open(name);
        zdd->dumpDot(file);
        file.close();
    }

    int get_remove()
    {
        return remove_nodes;
    }

    void init_table_farkas()
    {
        tdzdd::NodeTableHandler<2>& node_handler_farkas = zdd->getDiagram();
        farkas_table.init(njobs + 1);
        size_t const m = node_handler_farkas.privateEntity()[0].size();
        farkas_table[0].resize(m);

        for (size_t i = 0; i < m; ++i) {
            farkas_table[0][i].init_terminal_node(i);
        }

        for (int i = 1; i <= njobs; ++i) {
            tdzdd::MyVector<tdzdd::Node<2>> const& node =
                                             node_handler_farkas.privateEntity()[i];
            size_t const mm = node.size();
            farkas_table[i].resize(mm);

            for (size_t j = 0; j < mm; ++j) {
                farkas_table[i][j].init_node();
            }
        }
    }

    void evaluate_nodes(double *pi, int UB, double LB, int nmachines,
                        double reduced_cost)
    {
        tdzdd::NodeId&              root = zdd->root();

        /** Calculate the distance from  the origin to the given node */
        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();
            int          layer = nlayers - i;
            job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                              interval_list, layer);
            Job *job = tmp_pair->j;

            for (size_t j = 0; j < m; j++) {
                for (my_iterator<double> it = zdd_table[i][j].list.begin();
                        it != zdd_table[i][j].list.end(); it++) {
                    if (i == root.row()) {
                        (*it)->dist = 0;
                    }

                    if ((*it)->y->dist < (*it)->dist + pi[job->job] - value_Fj((
                                *it)->weight + job->processingime, job)) {
                        (*it)->y->dist = (*it)->dist + pi[job->job] - value_Fj((*it)->weight +
                                         job->processingime, job);
                    }

                    if ((*it)->n->dist < (*it)->dist) {
                        (*it)->n->dist  = (*it)->dist;
                    }
                }
            }
        }

        /** check for each node the Lagrangian dual */
        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();

            for (size_t j = 0; j < m; j++) {
                for (my_iterator<double> it = zdd_table[i][j].list.begin();
                        it != zdd_table[i][j].list.end(); it++) {
                    if (ceil(LB) - (double)(nmachines - 1)*reduced_cost - ((*it)->dist +
                            (*it)->b) > UB - 1 && ((*it)->calc)) {
                        (*it)->calc = false;
                        remove_edges++;
                    }

                    if(ceil(LB) - (double)(nmachines - 1)*reduced_cost - ((*it)->dist +
                            (*it)->c) > UB - 1 && ((*it)->calc0)) {
                        (*it)->calc0 = false;
                        remove_edges++;
                    }

                    if((*it)->calc0 == false && (*it)->calc == false) {
                        remove_nodes++;
                    }
                }
            }
        }


        printf("removed edges = %d, removed nodes = %d\n", remove_edges, remove_nodes);
    }

    inline bool evaluate_edge(tdzdd::NodeId *node)
    {
        return (node->row() > 0 || (node->row() == 0 && node->col() > 0));
    }

    void init_zdd_table()
    {
        tdzdd::NodeTableHandler<2>& handler = zdd->getDiagram();
        tdzdd::NodeId&              root = zdd->root();
        zdd_table.init(root.row() + 1);

        /** init table */
        for (int i = root.row(); i >= 0; i--) {
            tdzdd::MyVector<tdzdd::Node<2>> const& node = handler.privateEntity()[i];
            size_t const m = node.size();
            zdd_table[i].resize(m);
        }

        /** init root */
        zdd_table[root.row()][root.col()].add_weight(0, 0);

        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();
            int          layer = nlayers - i;
            job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                              interval_list, layer);
            Job *job = tmp_pair->j;

            for (size_t j = 0; j < m; j++) {
                tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
                tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

                for (my_iterator<double> it = zdd_table[i][j].list.begin();
                        it != zdd_table[i][j].list.end(); it++) {
                    (*it)->y = zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight((
                                   *it)->weight + job->processingime, nlayers - cur_node_1.row());

                    if ((cur_node_1.row() > 0) || (cur_node_1.row() == 0 && cur_node_1.col() == 1U)) {
                        edges.push_back(edge<double>(value_Fj((*it)->weight + job->processingime, job),
                                                     job, *it, (*it)->y));
                    }

                    (*it)->n = zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight((*it)->weight,
                               nlayers - cur_node_0.row());

                    if ((cur_node_0.row() > 0) || (cur_node_0.row() == 0 && cur_node_0.col() == 1U)) {
                        edges.push_back(edge<double>(0.0, NULL, *it, (*it)->n));
                    }
                }
            }
        }

        printf("number of edges %ld\n", edges.size());
        /** init terminal nodes */
        size_t const mm = handler.privateEntity()[0].size();

        for (size_t j = 0; j < mm; j++) {
            zdd_table[0][j].init_terminal_node(j);
        }
    }

    void init_zdd_one_conflict(int v1, int v2, int same)
    {
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
                                  int  ecount_differ)
    {
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

    void free_zdd_solver(int ecount_same, int ecount_differ)
    {
        if (ecount_same + ecount_differ > 0) {
            delete zdd;
        }

        zdd_table.init();
        farkas_table.init();
    }


    class Optimal_Solution<double>
        solve_weight_zdd_double(double *pi)
    {
        return zdd->evaluate_weight(WeightZDDdouble(pi, interval_list, njobs), zdd_table);
    }

    class Optimal_Solution<double>
        solve_farkas_double(double *pi)
    {
        return Optimal_Solution<double>();
    }

    void iterate_zdd()
    {
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
