#include <PricerSolver.hpp>
#include <set>
#include <vector>

/**
 * PricerSolverBase default COnstructor
 */
PricerSolverBase::PricerSolverBase(GPtrArray *_jobs):
    jobs(_jobs), njobs(_jobs->len), ordered_jobs(nullptr), nlayers(0) {
    dd = nullptr;
    zdd = nullptr;
    nb_nodes_bdd = 0;
    nb_nodes_zdd = 0;
}


PricerSolverBase::PricerSolverBase(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    jobs(_jobs), njobs(_jobs->len), ordered_jobs(_ordered_jobs),
    nlayers(ordered_jobs->len) {
    /**
     * Construction of decision diagram
     */
    PricerConstruct ps(ordered_jobs);
    dd = new tdzdd::DdStructure<2>(ps);
    nb_nodes_bdd = dd->size();

    /**
     * Construction of ZDD
     */
    zdd = new tdzdd::DdStructure<2>;
    *zdd = *dd;
    zdd->zddReduce();
    nb_nodes_zdd = zdd->size();
}

PricerSolverBase::~PricerSolverBase() {
    if (dd) {
        delete dd;
    }

    if (zdd) {
        delete zdd;
    }
}

/**
 * Some getters
 */
void PricerSolverBase::IterateZdd() {
    tdzdd::DdStructure<2>::const_iterator it = zdd->begin();

    for (; it != zdd->end(); ++it) {
        std::set<int>::const_iterator i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << nlayers - *i << " ";
        }

        std::cout << std::endl;
    }
}

void PricerSolverBase::PrintNumberPaths() {
    cout << "Number of paths: " << zdd->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

void PricerSolverBase::create_dot_zdd(const char *name) {
    std::ofstream file;
    file.open(name);
    zdd->dumpDot(file);
    file.close();
}

void PricerSolverBase::print_number_nodes_edges() {
    printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
           nb_removed_nodes);
}

int PricerSolverBase::get_remove() {
    return nb_removed_nodes;
}

size_t PricerSolverBase::get_datasize() {
    return dd->size();
}

size_t PricerSolverBase::get_numberrows_zdd() {
    return zdd->root().row();
}

double PricerSolverBase::get_cost_edge(int idx) {
    return 0.0;
}

/**
 * Reduced cost Fixing
 */
void PricerSolverBase::evaluate_nodes(double *pi, int UB, double LB,
                                      int nmachines, double reduced_cost) {
    return;
}

void PricerSolverBase::calculate_new_ordered_jobs() {
    return;
}

void PricerSolverBase::calculate_edges(scheduleset *set) {
    return;
}

/**
 * PricerSolverBdd constructor
 */
PricerSolverBdd::PricerSolverBdd(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBase(_jobs, _ordered_jobs) {
    InitTable();
}

void PricerSolverBdd::InitTable() {
    tdzdd::NodeTableHandler<2> &handler = dd->getDiagram();
    table.init(dd->topLevel() + 1);

    /** init table */
    for (int i = zdd->topLevel(); i >= 0; i--) {
        tdzdd::MyVector<tdzdd::Node<2>> const &layer = handler.privateEntity()[i];
        size_t const m = layer.size();
        table[i].resize(m);

        for (auto &it : table[i]) {
            if (i != 0) {
                int          layer = nlayers - i;
                job_interval_pair *tmp_pair = reinterpret_cast<job_interval_pair *>
                                              (g_ptr_array_index(
                                                   ordered_jobs, layer));
                it.set_job(tmp_pair->j);
            } else {
                it.set_job(nullptr, true);
            }
        }
    }

    /** init root */
    table[zdd->topLevel()][0].InitNode(0, true);

    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = table[i].size();

        for (size_t j = 0; j < m; j++) {
            Node<double> &n = table[i][j];
            Job *job = n.GetJob();

            tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
            tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

            n.child[1] = table[cur_node_1.row()][cur_node_1.col()].InitNode(
                             n.GetWeight() + job->processingime);

            n.child[0] = table[cur_node_0.row()][cur_node_0.col()].InitNode(n.GetWeight());

            // if ((cur_node_1.row() > 0) || (cur_node_1.row() == 0 &&
            //                                cur_node_1.col() == 1U)) {
            //     edges.push_back(std::make_shared<edge<double>>(value_Fj(it->weight +
            //                     job->processingime, job), job, it, it->y));
            //     it->out_edge.push_back(edges.back());
            //     it->y->in_edge.push_back(edges.back());
            // }

            // if ((cur_node_0.row() > 0) || (cur_node_0.row() == 0 &&
            //                                cur_node_0.col() == 1U)) {
            //     edges.push_back(make_shared<edge<double>>(0.0, nullptr, it, it->n));
            //     it->out_edge.push_back(edges.back());
            //     it->n->in_edge.push_back(edges.back());
            // }
        }
    }
}

PricerSolverBddSimple::PricerSolverBddSimple(GPtrArray *_jobs,
        GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _ordered_jobs) {
    evaluator = ForwardBddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return dd->evaluate_forward(&evaluator, table);
}

PricerSolverBddCycle::PricerSolverBddCycle(GPtrArray *_jobs,
        GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _ordered_jobs) {
    evaluator = ForwardBddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return dd->evaluate_forward(&evaluator, table);
}

/**
 * PricerSolverZdd Constructor
 */
PricerSolverZdd::PricerSolverZdd(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBase(_jobs, _ordered_jobs) {
    InitTable();
}

void PricerSolverZdd::InitTable() {
    tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
    table.init(zdd->topLevel() + 1);

    /** Init table */
    for (int i = zdd->topLevel(); i >= 0; i--) {
        tdzdd::MyVector<tdzdd::Node<2>> const &layer = handler.privateEntity()[i];
        size_t const m = layer.size();
        table[i].resize(m);

        for (auto &it : table[i]) {
            if (i != 0) {
                int          layer = nlayers - i;
                job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                                  ordered_jobs, layer);
                it.set_job(tmp_pair->j);
            } else {
                it.set_job(nullptr);
            }
        }
    }

    /**
     * Init root
     */
    table[zdd->topLevel()][0].add_weight(0, 0);

    /**
     * Construct in top-down fashion all the nodes in the diagram
     */
    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = table[i].size();

        for (size_t j = 0; j < m; j++) {
            Job *job = table[i][j].get_job();
            tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
            tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

            for (auto &it : table[i][j].list) {
                it->y = table[cur_node_1.row()][cur_node_1.col()].add_weight((
                            it)->GetWeight() + job->processingime, nlayers - cur_node_1.row());

                it->n = table[cur_node_0.row()][cur_node_0.col()].add_weight(it->GetWeight(),
                        nlayers - cur_node_0.row());
            }
        }
    }
}


PricerSolverCycle::PricerSolverCycle(GPtrArray *_jobs,
                                     GPtrArray *_ordered_jobs) :
    PricerSolverZdd(_jobs, _ordered_jobs) {
    evaluator = ForwardZddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return zdd->evaluate_forward(&evaluator, table);
}

PricerSolverZddSimple::PricerSolverZddSimple(GPtrArray *_jobs,
        GPtrArray *_ordered_jobs) :
    PricerSolverZdd(_jobs, _ordered_jobs) {
    evaluator = ForwardZddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverZddSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return zdd->evaluate_forward(&evaluator, table);
}

PricerSolverSimpleDp::PricerSolverSimpleDp(GPtrArray *_jobs, int _Hmax):
    PricerSolverBase(_jobs), Hmax(_Hmax) {}

void PricerSolverSimpleDp::InitTable() {
}

Optimal_Solution<double> PricerSolverSimpleDp::pricing_algorithm(double *_pi) {
    Optimal_Solution<double> opt_sol;
    opt_sol.cost = 0;
    double *F;
    Job **A;
    int t_min = 0;
    F = new double[Hmax + 1];
    A = new Job* [Hmax + 1];
    std::vector<Job *> v;


    /** Initialisation */
    F[0] = _pi[njobs];
    A[0] = nullptr;

    for (int t = 1; t < Hmax + 1; t++) {
        F[t] = -DBL_MAX / 2;
        A[t] = nullptr;
    }


    /** Recursion */
    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < njobs + 1; i++) {
            int j = i - 1;
            Job *job = reinterpret_cast<Job *>(g_ptr_array_index(jobs, j));

            if (t >=  job->processingime) {
                if (F[t - job->processingime] - (double) value_Fj(t,
                        job) + _pi[job->job] >= F[t]) {
                    F[t] = F[t - job->processingime] - value_Fj(t, job) + _pi[job->job];
                    A[t] = job;
                }
            }
        }
    }

    /** Find optimal solution */
    opt_sol.obj = -DBL_MAX;

    for (int i =  0; i < Hmax + 1; i++) {
        if (F[i] > opt_sol.obj) {
            opt_sol.C_max = i;
            opt_sol.obj = F[i];
        }
    }

    t_min = opt_sol.C_max;

    /** Construct the solution */
    while (A[t_min] != nullptr) {
        Job *job = A[t_min];
        v.push_back(A[t_min]);
        opt_sol.cost += value_Fj(t_min, A[t_min]);
        t_min -= job->processingime;
    }

    std::vector<Job *>::reverse_iterator it = v.rbegin();

    for (; it != v.rend(); ++it) {
        g_ptr_array_add(opt_sol.jobs, *it) ;
    }

    /** Free the memory */
    delete[] A;
    delete[] F;
    return opt_sol;
}

/**
 * Default Constructor
 */

// /** Copy constructor */


// /** Move constructor */

// /** Move assignment operator */

// /** Destructor */

// void PricerSolver::build_mip(double *x_e) {
//     try {
//         printf("Building Mip model for the extented formulation:\n");
//         model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
//         model->set(GRB_IntParam_Threads, 1);
//         model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//         model->set(GRB_IntParam_Presolve, 2);
//         model->set(GRB_IntParam_VarBranch, 3);

//         /** adding variables */
//         for (auto &i : edges) {
//             if (i->job) {
//                 if (i->out->calc) {
//                     i->v = model->addVar(0.0, 1.0, i->cost, GRB_BINARY);
//                 } else {
//                     i->v = model->addVar(0.0, 0.0, i->cost, GRB_CONTINUOUS);
//                 }
//             } else {
//                 i->v = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
//             }
//         }

//         model->update();
//         GRBLinExpr obj = 0;

//         for (auto &i : edges) {
//             if (i->job && i->out->calc) {
//                 obj +=  i->cost * i->v;
//             }
//         }

//         model->addConstr(obj, GRB_LESS_EQUAL, (double) ub - 1.0);
//         model->update();
//         /** assign the jobs to some path */
//         printf("Adding assignment constraints:\n");

//         GRBLinExpr *assignment = new GRBLinExpr[njobs];

//         for (unsigned i = 0; i < jobs->len; ++i) {
//             assignment[i] = 0;
//         }

//         for (const auto &i : edges) {
//             if (i->job) {
//                 assignment[i->job->job] += i->v;
//             }
//         }

//         for (unsigned i = 0; i < jobs->len; ++i) {
//             model->addConstr(assignment[i], GRB_EQUAL, 1.0);
//         }

//         delete[] assignment;

//         model->update();

//         tdzdd::NodeId              &root = zdd->root();
//         /** Calculate the distance from  the origin to the given node */
//         printf("Adding flow constraints:\n");

//         for (int i = root.row() - 1; i > 0; i--) {
//             size_t const m = zdd_table[i].size();

//             for (size_t j = 0; j < m; j++) {
//                 for (auto &it : zdd_table[i][j].list) {
//                     GRBLinExpr expr = 0;

//                     for (const auto &v : it->out_edge) {
//                         auto p = v.lock();
//                         expr -= p->v;
//                     }

//                     for (const auto &v : it->in_edge) {
//                         auto p = v.lock();
//                         expr += p->v;
//                     }

//                     model->addConstr(expr, GRB_EQUAL, 0.0);
//                 }
//             }
//         }

//         printf("Adding convex constraint:\n");
//         GRBLinExpr expr = 0;

//         for (auto &it : zdd_table[0][1].list) {
//             for (const auto &v : it->in_edge) {
//                 auto p = v.lock();
//                 expr += p->v;
//             }
//         }

//         model->addConstr(expr, GRB_EQUAL, (double) nmachines);
//         model->update();
//         printf("Begin solving:\n");
//         model->optimize();
//     } catch (GRBException e) {
//         cout << "Error code = " << e.getErrorCode() << endl;
//         cout << e.getMessage() << endl;
//     } catch (...) {
//         cout << "Exception during optimization" << endl;
//     }
// }

// void PricerSolver::calculate_new_ordered_jobs() {
//     int          first_del = -1;
//     int          last_del = -1;
//     int it = 0;
//     tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();

//     /** remove the unnecessary edges of the zdd */
//     for (int i = zdd->topLevel(); i > 0; i--) {
//         size_t const m = zdd_table[i].size();
//         bool remove = true;

//         for (size_t j = 0; j < m; j++) {
//             for (const auto &iter : zdd_table[i][j].list) {
//                 if (iter->calc) {
//                     remove = false;
//                 } else {
//                     tdzdd::NodeId &cur_node_1 = handler.privateEntity().child(i, j, 1);
//                     cur_node_1 = 0;
//                 }
//             }
//         }

//         if (!remove) {
//             if (first_del != -1) {
//                 g_ptr_array_remove_range(ordered_jobs, first_del, last_del - first_del + 1);
//                 it = it - (last_del - first_del);
//                 first_del = last_del = -1;
//             } else {
//                 it++;
//             }
//         } else {
//             if (first_del == -1) {
//                 first_del = it;
//                 last_del = first_del;
//             } else {
//                 last_del++;
//             }

//             it++;
//         }
//     }

//     if (first_del != -1) {
//         g_ptr_array_remove_range(ordered_jobs, first_del, last_del - first_del + 1);
//     }

//     /** Construct the new zdd */
//     delete zdd;
//     nlayers = ordered_jobs->len;
//     PricerConstruct ps(ordered_jobs);
//     zdd = new tdzdd::DdStructure<2>(ps);
//     nb_nodes_bdd = zdd->size();
//     zdd->zddReduce();
//     int count = 0;

//     for (int i = njobs  - 1; i >= 0 && count <  8; --i) {
//         Job *tmp_j = (Job *) g_ptr_array_index(jobs, i);

//         if (tmp_j->nb_layers == 1) {
//             zdd->zddSubset(scheduling(tmp_j, ordered_jobs, 2));
//             zdd->zddReduce();
//             count++;
//         }
//     }

//     nb_nodes_zdd = zdd->size();
//     printf("The new number of layers = %u\n", ordered_jobs->len);
//     printf("The new size of BDD = %lu and size ZDD= %lu\n", nb_nodes_bdd,
//            nb_nodes_zdd);

//     /** calculate the new table */
//     zdd_table.init();
//     edges.clear();
//     init_zdd_table();
// }

// void PricerSolver::iterate_zdd() {
//     tdzdd::DdStructure<2>::const_iterator it = zdd->begin();

//     for (; it != zdd->end(); ++it) {
//         std::set<int>::const_iterator i = (*it).begin();

//         for (; i != (*it).end(); ++i) {
//             std::cout << nlayers - *i << " ";
//         }

//         std::cout << std::endl;
//     }
// }


// void PricerSolver::evaluate_nodes(double *pi, int UB, double LB, int nmachines,
//                                   double reduced_cost) {
//     double value;

//     /** Calculate the distance from  the origin to the given node */
//     for (int i = zdd->topLevel(); i > 0; i--) {
//         size_t const m = zdd_table[i].size();
//         int          layer = nlayers - i;
//         job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
//                                           ordered_jobs, layer);
//         Job *job = tmp_pair->j;

//         for (size_t j = 0; j < m; j++) {
//             for (auto &it : zdd_table[i][j].list) {
//                 if (i == zdd->topLevel()) {
//                     it->dist = 0;
//                 }

//                 value = pi[job->job] - value_Fj(it->weight + job->processingime, job);

//                 if (it->y->dist < it->dist + value) {
//                     it->y->dist = it->dist + value;
//                 }

//                 if (it->n->dist < it->dist) {
//                     it->n->dist  = it->dist;
//                 }
//             }
//         }
//     }

//     /** check for each node the Lagrangian dual */
//     for (int i = zdd->topLevel(); i > 0; i--) {
//         size_t const m = zdd_table[i].size();

//         for (size_t j = 0; j < m; j++) {
//             for (auto &it : zdd_table[i][j].list) {
//                 if (LB - (double)(nmachines - 1)*reduced_cost - (it->dist + it->b) > UB - 1 +
//                         0.0001 && (it->calc)) {
//                     it->calc = false;
//                     nb_removed_edges++;
//                 }

//                 if (LB - (double)(nmachines - 1)*reduced_cost - (it->dist + it->c) > UB - 1
//                         + 0.0001 && (it->calc0)) {
//                     it->calc0 = false;
//                     nb_removed_edges++;
//                 }

//                 if (it->calc0 == false && it->calc == false && it->remove_node == false) {
//                     nb_removed_nodes++;
//                     it->remove_node = true;
//                 }
//             }
//         }
//     }
// }

// Optimal_Solution<double> PricerSolver::dynamic_programming_ti(double *pi) {
//     Optimal_Solution<double> opt_sol;
//     opt_sol.cost = 0;
//     double *F;
//     Job **A;
//     int t_min = 0;
//     F = new double [Hmax + 1];
//     A = new Job* [Hmax + 1];
//     std::vector<Job *> v;


//     /** Initialisation */
//     F[0] = pi[njobs];
//     A[0] = nullptr;

//     for (int t = 1; t < Hmax + 1; t++) {
//         F[t] = -DBL_MAX / 2;
//         A[t] = nullptr;
//     }


//     /** Recursion */
//     for (int t = 0; t < Hmax + 1; t++) {
//         for (int i = 1; i < njobs + 1; i++) {
//             int j = i - 1;
//             Job *job = (Job *) g_ptr_array_index(jobs, j);

//             if (t >=  job->processingime) {
//                 if (F[t - job->processingime] - (double) value_Fj(t,
//                         job) + pi[job->job] >= F[t]) {
//                     F[t] = F[t - job->processingime] - value_Fj(t, job) + pi[job->job];
//                     A[t] = job;
//                 }
//             }
//         }
//     }

//     /** Find optimal solution */
//     opt_sol.obj = -DBL_MAX;

//     for (int i =  0; i < Hmax + 1; i++) {

//         if (F[i] > opt_sol.obj) {
//             opt_sol.C_max = i;
//             opt_sol.obj = F[i];
//         }
//     }

//     t_min = opt_sol.C_max;

//     /** Construct the solution */
//     while (A[t_min] != nullptr) {
//         Job *job = A[t_min];
//         v.push_back(A[t_min]);
//         opt_sol.cost += value_Fj(t_min, A[t_min]) ;
//         t_min -= job->processingime;
//     }

//     std::vector<Job *>::reverse_iterator it = v.rbegin();

//     for (; it != v.rend(); ++it) {
//         g_ptr_array_add(opt_sol.jobs, *it) ;
//     }

//     /** Free the memory */
//     delete[] A;
//     delete[] F;
//     return opt_sol;
// }