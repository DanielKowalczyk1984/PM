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
private:
    int nmachines;
    int ub;
    int njobs;
    int nlayers;
    GPtrArray *jobs;
    GPtrArray *ordered_jobs;
    tdzdd::DdStructure<2>                    *zdd;
    tdzdd::DataTable<PricerWeightZDD<double>> zdd_table;
    tdzdd::DataTable<PricerFarkasZDD<double>> farkas_table;
    std::vector<shared_ptr<edge<double>>> edges;
    int nb_removed_edges;
    int nb_removed_nodes;
    size_t nb_nodes_bdd;
    size_t nb_nodes_zdd;
    GRBEnv *env;
    GRBModel *model;


public:

    /** Default Constructor */
    PricerSolver(GPtrArray *_jobs, GPtrArray *_ordered_jobs, int _nmachines,
                 int _ub): nmachines(_nmachines), ub(_ub), njobs((int)_jobs->len),
        nlayers((int)_ordered_jobs->len), jobs(_jobs), ordered_jobs(_ordered_jobs)
    {
        /** Initialize some variables */
        nb_removed_edges = 0;
        nb_removed_nodes = 0;
        /** Construct the ZDD and the tables associated to it */
        PricerConstruct ps(ordered_jobs);
        zdd = new tdzdd::DdStructure<2>(ps);
        nb_nodes_bdd = zdd->size();
        zdd->zddReduce();
        nb_nodes_zdd = zdd->size();
        init_zdd_table();
        printf("size BDD = %lu, size ZDD= %lu\n", nb_nodes_bdd, nb_nodes_zdd);
        edges.reserve(nb_nodes_bdd);
        /** Initialize Gurobi Model */
        env = new GRBEnv();
        model = new GRBModel(*env);
    };

    /** Copy constructor */
    PricerSolver(const PricerSolver& other): nmachines(other.nmachines), ub(other.ub),
        njobs(other.njobs), nlayers(other.nlayers), jobs(other.jobs),
        ordered_jobs(other.ordered_jobs), zdd(new tdzdd::DdStructure<2>),
        zdd_table(other.zdd_table), farkas_table(other.farkas_table), edges(other.edges),
        nb_removed_edges(other.nb_removed_edges),
        nb_removed_nodes(other.nb_removed_nodes), nb_nodes_bdd(other.nb_nodes_bdd),
        nb_nodes_zdd(other.nb_nodes_zdd), env(new GRBEnv()), model(new GRBModel(*env))
    {
        zdd = other.zdd;
        env = other.env;
        model = other.model;
    }

    /** Move Constructor */
    PricerSolver(PricerSolver&& other) noexcept: nmachines(other.nmachines),
        ub(other.ub), njobs(other.njobs), nlayers(other.nlayers), jobs(other.jobs),
        ordered_jobs(other.ordered_jobs), zdd(other.zdd), zdd_table(other.zdd_table),
        farkas_table(other.farkas_table), edges(other.edges),
        nb_removed_edges(other.nb_removed_edges),
        nb_removed_nodes(other.nb_removed_nodes), nb_nodes_bdd(other.nb_nodes_bdd),
        nb_nodes_zdd(other.nb_nodes_zdd), env(other.env), model(other.model)
    {
        other.zdd = nullptr;
        other.env = nullptr;
        other.model = nullptr;
    }

    /** Move Constructor */
    PricerSolver& operator=(const PricerSolver& other)
    {
        PricerSolver tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    /** Move assignment operator */
    PricerSolver& operator=(PricerSolver&& other) noexcept
    {
        nmachines = other.nmachines;
        ub = other.ub;
        njobs = other.njobs;
        nlayers = other.nlayers;
        jobs = other.jobs;
        ordered_jobs = other.ordered_jobs;
        delete zdd;
        zdd = other.zdd;
        other.zdd = nullptr;
        zdd_table = other.zdd_table;
        farkas_table = other.farkas_table;
        edges = other.edges;
        nb_removed_edges = other.nb_removed_edges;
        nb_removed_nodes = other.nb_removed_nodes;
        nb_nodes_bdd = other.nb_nodes_bdd;
        nb_nodes_zdd = other.nb_nodes_zdd;
        /** Delete memory of Gurobi */
        delete model;
        delete env;
        env = other.env;
        model = other.model;
        other.env = nullptr;
        other.model = nullptr;
        return *this;
    }

    /** Destructor */
    ~PricerSolver() noexcept
    {
        delete zdd;
        delete model;
        delete env;
    }

    void build_mip()
    {
        try {
            printf("Building Mip model for the extented formulation:\n");
            model->set(GRB_IntParam_Method, GRB_METHOD_DUAL);
            model->set(GRB_IntParam_Threads, 1);
            model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model->set(GRB_IntParam_Presolve, 2);
            model->set(GRB_IntParam_VarBranch, 3);

            /** adding variables */
            for (auto& i : edges) {
                if (i->job) {
                    if (i->out->calc) {
                        i->v = model->addVar(0.0, 1.0, i->cost, GRB_BINARY);
                    } else {
                        i->v = model->addVar(0.0, 0.0, i->cost, GRB_CONTINUOUS);
                    }
                } else {
                    i->v = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
                }
            }

            model->update();
            GRBLinExpr obj = 0;

            for (auto& i : edges) {
                if (i->job && i->out->calc) {
                    obj +=  i->cost * i->v;
                }
            }

            model->addConstr(obj, GRB_LESS_EQUAL, (double) ub - 1.0);
            model->update();
            /** assign the jobs to some path */
            printf("Adding assignment constraints:\n");

            GRBLinExpr *assignment = new GRBLinExpr[njobs];

            for (unsigned i = 0; i < jobs->len; ++i) {
                assignment[i] = 0;
            }

            for (const auto& i : edges) {
                if (i->job) {
                    assignment[i->job->job] += i->v;
                }
            }

            for (unsigned i = 0; i < jobs->len; ++i) {
                model->addConstr(assignment[i], GRB_EQUAL, 1.0);
            }

            delete[] assignment;

            model->update();

            tdzdd::NodeId&              root = zdd->root();
            /** Calculate the distance from  the origin to the given node */
            printf("Adding flow constraints:\n");
            for (int i = root.row() - 1; i > 0; i--) {
                size_t const m = zdd_table[i].size();
                for (size_t j = 0; j < m; j++) {
                    for (auto& it : zdd_table[i][j].list) {
                        GRBLinExpr expr = 0;

                        for (const auto& v : it->out_edge) {
                                expr -= v->v;
                        }

                        for(const auto& v : it->in_edge) {
                            expr += v->v;
                        }

                        model->addConstr(expr, GRB_EQUAL, 0.0);
                    }
                }
            }

            printf("Adding convex constraint:\n");
            GRBLinExpr expr = 0;

            for (auto& it : zdd_table[0][1].list) {
                for (const auto& v : it->in_edge) {
                    expr += v->v;
                }
            }

            model->addConstr(expr, GRB_EQUAL, (double) nmachines);
            model->update();
            printf("Begin solving:\n");
            model->optimize();
        } catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Exception during optimization" << endl;
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
        return nb_removed_nodes;
    }

    size_t get_datasize()
    {
        return zdd->size();
    }

    size_t get_numberrows_zdd()
    {
        return zdd->root().row();
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

    void calculate_new_ordered_jobs()
    {
        tdzdd::NodeId&              root = zdd->root();
        int          first_del = -1;
        int          last_del = -1;
        int it = 0;

        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();
            bool remove = true;

            for (size_t j = 0; j < m && remove; j++) {
                for (const auto& iter : zdd_table[i][j].list) {
                    if (iter->calc) {
                        remove = false;
                    }
                }
            }

            if (!remove) {
                if (first_del != -1) {
                    g_ptr_array_remove_range(ordered_jobs, first_del, last_del - first_del + 1);
                    it = it - (last_del - first_del);
                    first_del = last_del = -1;
                } else {
                    it++;
                }
            } else {
                if (first_del == -1) {
                    first_del = it;
                    last_del = first_del;
                } else {
                    last_del++;
                }

                it++;
            }
        }

        if (first_del != -1) {
            g_ptr_array_remove_range(ordered_jobs, first_del, last_del - first_del + 1);
        }

        printf("The new number of layers = %u\n", ordered_jobs->len);
        delete zdd;
        nlayers = ordered_jobs->len;
        PricerConstruct ps(ordered_jobs);
        zdd = new tdzdd::DdStructure<2>(ps);
        nb_nodes_bdd = zdd->size();
        zdd->zddReduce();
        nb_nodes_zdd = zdd->size();
        printf("size BDD = %lu, size ZDD= %lu\n", nb_nodes_bdd, nb_nodes_zdd);
        zdd_table.init();
        edges.clear();
        init_zdd_table();
    }

    void evaluate_nodes(double *pi, int UB, double LB, int nmachines,
                        double reduced_cost)
    {
        tdzdd::NodeId&              root = zdd->root();
        double value;

        /** Calculate the distance from  the origin to the given node */
        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();
            int          layer = nlayers - i;
            job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                              ordered_jobs, layer);
            Job *job = tmp_pair->j;

            for (size_t j = 0; j < m; j++) {
                for (auto& it : zdd_table[i][j].list) {
                    if (i == root.row()) {
                        it->dist = 0;
                    }

                    value = pi[job->job] - value_Fj(it->weight +job->processingime, job);
                    if (it->y->dist < it->dist + value ) {
                        it->y->dist = it->dist + value;
                    }

                    if (it->n->dist < it->dist) {
                        it->n->dist  = it->dist;
                    }
                }
            }
        }

        /** check for each node the Lagrangian dual */
        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();

            for (size_t j = 0; j < m; j++) {
                for (auto& it : zdd_table[i][j].list) {
                    if (LB - (double)(nmachines - 1)*reduced_cost - (it->dist + it->b) > UB - 1 +
                            0.0001 && (it->calc)) {
                        it->calc = false;
                        nb_removed_edges++;
                    }

                    if (LB - (double)(nmachines - 1)*reduced_cost - (it->dist + it->c) > UB - 1
                            && (it->calc0)) {
                        it->calc0 = false;
                        nb_removed_edges++;
                    }

                    if (it->calc0 == false && it->calc == false && it->remove_node == false) {
                        nb_removed_nodes++;
                        it->remove_node = true;
                    }
                }
            }
        }
    }

    void print_number_nodes_edges()
    {
        printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
               nb_removed_nodes);
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
        zdd_table[root.row()][root.col()].add_weight(0, 0, njobs);

        for (int i = root.row(); i > 0; i--) {
            size_t const m = zdd_table[i].size();
            int          layer = nlayers - i;
            job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                              ordered_jobs, layer);
            Job *job = tmp_pair->j;

            for (size_t j = 0; j < m; j++) {
                tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
                tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

                for (auto& it : zdd_table[i][j].list) {
                    it->y = zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight((
                                it)->weight + job->processingime, nlayers - cur_node_1.row(), njobs);

                    if ((cur_node_1.row() > 0) || (cur_node_1.row() == 0 && cur_node_1.col() == 1U)) {
                        edges.push_back(std::make_shared<edge<double>>(value_Fj(it->weight +
                                        job->processingime, job), job, it, it->y));
                    }
                    it->out_edge.push_back(edges.back());
                    it->y->in_edge.push_back(edges.back());

                    it->n = zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight(it->weight,
                            nlayers - cur_node_0.row(), njobs);

                    if ((cur_node_0.row() > 0) || (cur_node_0.row() == 0 && cur_node_0.col() == 1U)) {
                        edges.push_back(make_shared<edge<double>>(0.0, nullptr, it, it->n));
                    }
                    it->out_edge.push_back(edges.back());
                    it->n->in_edge.push_back(edges.back());
                }
            }
        }

        /** init terminal nodes */
        size_t const mm = handler.privateEntity()[0].size();

        for (size_t j = 0; j < mm; j++) {
            zdd_table[0][j].init_terminal_node(j, njobs);
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


    class Optimal_Solution<double> solve_weight_zdd_double(double *pi)
    {
        return zdd->evaluate_weight(WeightZDDdouble(pi, ordered_jobs, njobs), zdd_table);
    }

    class Optimal_Solution<double> solve_farkas_double(double *pi)
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
