#include <PricerSolver.hpp>

/**
 * Default Constructor
 */
PricerSolver::PricerSolver(GPtrArray *_jobs, GPtrArray *_ordered_jobs,
                           int _nmachines,
                           int _ub, int _Hmax): nmachines(_nmachines), ub(_ub), Hmax(_Hmax),
    njobs((int)_jobs->len),
    nlayers((int)_ordered_jobs->len), jobs(_jobs), ordered_jobs(_ordered_jobs) {
    /** Initialize some variables */
    nb_removed_edges = 0;
    nb_removed_nodes = 0;
    /** Construct the ZDD and the tables associated to it */
    PricerConstruct ps(ordered_jobs);
    dd = new tdzdd::DdStructure<2>(ps);
    nb_nodes_bdd = dd->size();
    zdd = new tdzdd::DdStructure<2>;
    *zdd = *dd;
    zdd->zddReduce();
    nb_nodes_zdd = zdd->size();
    init_zdd_table();
    init_bdd_table();
    init_zdd_duration_table();
    printf("The size of BDD = %lu and size ZDD= %lu\n", nb_nodes_bdd, nb_nodes_zdd);
    edges.reserve(2 * nb_nodes_bdd);
    /** Initialize Gurobi Model */
    env = new GRBEnv();
    model = new GRBModel(*env);
}

/** Copy constructor */
PricerSolver::PricerSolver(const PricerSolver &other): nmachines(
        other.nmachines), ub(other.ub),
    njobs(other.njobs), nlayers(other.nlayers), jobs(other.jobs),
    ordered_jobs(other.ordered_jobs), zdd(new tdzdd::DdStructure<2>),
    dd(new tdzdd::DdStructure<2>),
    zdd_table(other.zdd_table), farkas_table(other.farkas_table),
    edges(other.edges),
    nb_removed_edges(other.nb_removed_edges),
    nb_removed_nodes(other.nb_removed_nodes), nb_nodes_bdd(other.nb_nodes_bdd),
    nb_nodes_zdd(other.nb_nodes_zdd), env(new GRBEnv()), model(new GRBModel(*env)) {
    dd = other.dd;
    zdd = other.zdd;
    env = other.env;
    model = other.model;
}

/** Move constructor */
PricerSolver::PricerSolver(PricerSolver &&other) noexcept: nmachines(
        other.nmachines),
    ub(other.ub), njobs(other.njobs), nlayers(other.nlayers), jobs(other.jobs),
    ordered_jobs(other.ordered_jobs), zdd(other.zdd), dd(other.dd),
    zdd_table(other.zdd_table),
    farkas_table(other.farkas_table), edges(other.edges),
    nb_removed_edges(other.nb_removed_edges),
    nb_removed_nodes(other.nb_removed_nodes), nb_nodes_bdd(other.nb_nodes_bdd),
    nb_nodes_zdd(other.nb_nodes_zdd), env(other.env), model(other.model) {
    other.dd = nullptr;
    other.zdd = nullptr;
    other.env = nullptr;
    other.model = nullptr;
}

PricerSolver &PricerSolver::operator=(const PricerSolver &other) {
    PricerSolver tmp(other);
    *this = std::move(tmp);
    return *this;
}

/** Move assignment operator */
PricerSolver &PricerSolver::operator=(PricerSolver &&other) noexcept {
    if (this != &other) {
        nmachines = other.nmachines;
        ub = other.ub;
        njobs = other.njobs;
        nlayers = other.nlayers;
        jobs = other.jobs;
        ordered_jobs = other.ordered_jobs;
        delete zdd;
        zdd = other.zdd;
        delete dd;
        dd = other.dd;
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
    }

    return *this;
}

/** Destructor */
PricerSolver::~PricerSolver() noexcept {
    delete dd;
    delete zdd;
    delete model;
    delete env;
}

void PricerSolver::build_mip(double *x_e) {
    try {
        printf("Building Mip model for the extented formulation:\n");
        model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
        model->set(GRB_IntParam_Threads, 1);
        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model->set(GRB_IntParam_Presolve, 2);
        model->set(GRB_IntParam_VarBranch, 3);

        /** adding variables */
        for (auto &i : edges) {
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

        for (auto &i : edges) {
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

        for (const auto &i : edges) {
            if (i->job) {
                assignment[i->job->job] += i->v;
            }
        }

        for (unsigned i = 0; i < jobs->len; ++i) {
            model->addConstr(assignment[i], GRB_EQUAL, 1.0);
        }

        delete[] assignment;

        model->update();

        tdzdd::NodeId              &root = zdd->root();
        /** Calculate the distance from  the origin to the given node */
        printf("Adding flow constraints:\n");

        for (int i = root.row() - 1; i > 0; i--) {
            size_t const m = zdd_table[i].size();

            for (size_t j = 0; j < m; j++) {
                for (auto &it : zdd_table[i][j].list) {
                    GRBLinExpr expr = 0;

                    for (const auto &v : it->out_edge) {
                        auto p = v.lock();
                        expr -= p->v;
                    }

                    for (const auto &v : it->in_edge) {
                        auto p = v.lock();
                        expr += p->v;
                    }

                    model->addConstr(expr, GRB_EQUAL, 0.0);
                }
            }
        }

        printf("Adding convex constraint:\n");
        GRBLinExpr expr = 0;

        for (auto &it : zdd_table[0][1].list) {
            for (const auto &v : it->in_edge) {
                auto p = v.lock();
                expr += p->v;
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

void PricerSolver::calculate_new_ordered_jobs() {
    int          first_del = -1;
    int          last_del = -1;
    int it = 0;
    tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();

    /** remove the unnecessary edges of the zdd */
    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = zdd_table[i].size();
        bool remove = true;

        for (size_t j = 0; j < m; j++) {
            for (const auto &iter : zdd_table[i][j].list) {
                if (iter->calc) {
                    remove = false;
                } else {
                    tdzdd::NodeId &cur_node_1 = handler.privateEntity().child(i, j, 1);
                    cur_node_1 = 0;
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

    /** Construct the new zdd */
    delete zdd;
    nlayers = ordered_jobs->len;
    PricerConstruct ps(ordered_jobs);
    zdd = new tdzdd::DdStructure<2>(ps);
    nb_nodes_bdd = zdd->size();
    zdd->zddReduce();
    int count = 0;

    for (int i = njobs  - 1; i >= 0 && count <  8; --i) {
        Job *tmp_j = (Job *) g_ptr_array_index(jobs, i);

        if (tmp_j->nb_layers == 1) {
            zdd->zddSubset(scheduling(tmp_j, ordered_jobs, 2));
            zdd->zddReduce();
            count++;
        }
    }

    nb_nodes_zdd = zdd->size();
    printf("The new number of layers = %u\n", ordered_jobs->len);
    printf("The new size of BDD = %lu and size ZDD= %lu\n", nb_nodes_bdd,
           nb_nodes_zdd);

    /** calculate the new table */
    zdd_table.init();
    edges.clear();
    init_zdd_table();
}

void PricerSolver::iterate_zdd() {
    tdzdd::DdStructure<2>::const_iterator it = zdd->begin();

    for (; it != zdd->end(); ++it) {
        std::set<int>::const_iterator i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << nlayers - *i << " ";
        }

        std::cout << std::endl;
    }
}

void PricerSolver::print_number_paths() {
    cout << "Number of paths: " << zdd->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

void PricerSolver::print_number_nodes_edges() {
    printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
           nb_removed_nodes);
}

void PricerSolver::evaluate_nodes(double *pi, int UB, double LB, int nmachines,
                    double reduced_cost) {
    double value;

    /** Calculate the distance from  the origin to the given node */
    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = zdd_table[i].size();
        int          layer = nlayers - i;
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          ordered_jobs, layer);
        Job *job = tmp_pair->j;

        for (size_t j = 0; j < m; j++) {
            for (auto &it : zdd_table[i][j].list) {
                if (i == zdd->topLevel()) {
                    it->dist = 0;
                }

                value = pi[job->job] - value_Fj(it->weight + job->processingime, job);

                if (it->y->dist < it->dist + value) {
                    it->y->dist = it->dist + value;
                }

                if (it->n->dist < it->dist) {
                    it->n->dist  = it->dist;
                }
            }
        }
    }

    /** check for each node the Lagrangian dual */
    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = zdd_table[i].size();

        for (size_t j = 0; j < m; j++) {
            for (auto &it : zdd_table[i][j].list) {
                if (LB - (double)(nmachines - 1)*reduced_cost - (it->dist + it->b) > UB - 1 +
                        0.0001 && (it->calc)) {
                    it->calc = false;
                    nb_removed_edges++;
                }

                if (LB - (double)(nmachines - 1)*reduced_cost - (it->dist + it->c) > UB - 1
                        + 0.0001 && (it->calc0)) {
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

void PricerSolver::init_zdd_table() {
    tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
    zdd_table.init(zdd->topLevel() + 1);
    reset_total();

    /** init table */
    for (int i = zdd->topLevel(); i >= 0; i--) {
        tdzdd::MyVector<tdzdd::Node<2>> const &layer = handler.privateEntity()[i];
        size_t const m = layer.size();
        zdd_table[i].resize(m);
    }

    /** init root */
    zdd_table[zdd->topLevel()][0].add_weight(0, 0, njobs);

    /** init terminal nodes */
    size_t const mm = zdd_table[0].size();

    for (unsigned i = 0; i < mm; ++i) {
        zdd_table[0][i].add_terminal_node(i, nlayers, njobs);
    }

    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = zdd_table[i].size();

        int          layer = nlayers - i;
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          ordered_jobs, layer);
        Job *job = tmp_pair->j;

        for (size_t j = 0; j < m; j++) {
            tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
            tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

            for (auto &it : zdd_table[i][j].list) {
                it->y = zdd_table[cur_node_1.row()][cur_node_1.col()].add_weight((
                            it)->weight + job->processingime, nlayers - cur_node_1.row(), njobs);
                
                it->n = zdd_table[cur_node_0.row()][cur_node_0.col()].add_weight(it->weight,
                        nlayers - cur_node_0.row(), njobs);

                if ((cur_node_1.row() > 0) || (cur_node_1.row() == 0 &&
                                               cur_node_1.col() == 1U)) {
                    edges.push_back(std::make_shared<edge<double>>(value_Fj(it->weight +
                                    job->processingime, job), job, it, it->y));
                    it->out_edge.push_back(edges.back());
                    it->y->in_edge.push_back(edges.back());
                }

                if ((cur_node_0.row() > 0) || (cur_node_0.row() == 0 &&
                                               cur_node_0.col() == 1U)) {
                    edges.push_back(make_shared<edge<double>>(0.0, nullptr, it, it->n));
                    it->out_edge.push_back(edges.back());
                    it->n->in_edge.push_back(edges.back());
                }
            }
        }
    }
}

void PricerSolver::init_bdd_table() {
    tdzdd::NodeTableHandler<2> &handler = dd->getDiagram();
    tdzdd::NodeId &root = dd->root();
    dd_table.init(root.row() + 1);

    /** init table */
    for (int i = root.row(); i >= 0 ; i--) {
        tdzdd::MyVector<tdzdd::Node<2>> const &node = handler.privateEntity()[i];
        size_t const m = node.size();
        dd_table[i].resize(m);
    }

    /** init root */
    dd_table[root.row()][root.col()].init_node(0, njobs);

    for (size_t i = root.row(); i > 0 ; i--) {
        size_t const m = dd_table[i].size();
        int layer = nlayers - i;
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          ordered_jobs, layer);
        Job *job = tmp_pair->j;

        for (size_t j = 0; j < m; j++) {
            int sum_p = dd_table[i][j].sum_p;
            tdzdd::NodeId cur_node = handler.privateEntity().child(i, j, 0);

            if (cur_node.row() != 0) {
                dd_table[cur_node.row()][cur_node.col()].init_node(sum_p, njobs);
            }

            cur_node = handler.privateEntity().child(i, j, 1);

            if (cur_node.row() != 0) {
                dd_table[cur_node.row()][cur_node.col()].init_node(sum_p + job->processingime,
                        njobs);
            }
        }
    }

    /** init terminal nodes */
    size_t const mm = handler.privateEntity()[0].size();

    for (size_t j  = 0; j < mm; j++) {
        dd_table[0][j].init_terminal_node(j);
    }
}

void PricerSolver::init_table_farkas() {
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

void PricerSolver::init_zdd_duration_table() {
    tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
    zdd_duration_table.init(zdd->topLevel() + 1);

    /** init table */
    for (int i = zdd->topLevel(); i >= 0; i--) {
        tdzdd::MyVector<tdzdd::Node<2>> const &layer = handler.privateEntity()[i];
        size_t const m = layer.size();
        zdd_duration_table[i].resize(m);
    }

    /** init root */
    zdd_duration_table[zdd->topLevel()][0].add_weight(0, 0,true, false);

    for (int i = zdd->topLevel(); i > 0; i--) {
        size_t const m = zdd_duration_table[i].size();

        int          layer = nlayers - i;
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(ordered_jobs, layer);
        Job *job = tmp_pair->j;

        for (size_t j = 0; j < m; j++) {
            tdzdd::NodeId cur_node_0 = handler.privateEntity().child(i, j, 0);
            tdzdd::NodeId cur_node_1 = handler.privateEntity().child(i, j, 1);

            for (NodeDuration<double> &it : zdd_duration_table[i][j].list) {
                if(cur_node_1.row() != 0) {
                    it.y = zdd_duration_table[cur_node_1.row()][cur_node_1.col()].add_weight(it.GetWeight() + job->processingime, nlayers - cur_node_1.row());
                } else if (cur_node_1.col() == 1) {
                    it.y = zdd_duration_table[cur_node_1.row()][cur_node_1.col()].add_weight(it.GetWeight() + job->processingime, nlayers, false, true); 
                }

                if(cur_node_0.row() != 0) {
                    it.n = zdd_duration_table[cur_node_0.row()][cur_node_0.col()].add_weight(it.GetWeight(), nlayers - cur_node_0.row());
                } else if (cur_node_0.code() == 1) {
                    it.n = zdd_duration_table[cur_node_0.row()][cur_node_0.col()].add_weight(it.GetWeight(), nlayers, false, true);
                }
                it.SetJob(job);
            }
        }
    }
}

void PricerSolver::init_tables() {
    init_zdd_table();
    init_table_farkas();
    init_zdd_duration_table();
}


void PricerSolver::create_dot_zdd(const char *name) {
    std::ofstream file;
    file.open(name);
    zdd->dumpDot(file);
    file.close();
}

int PricerSolver::get_remove() {
    return nb_removed_nodes;
}

size_t PricerSolver::get_datasize() {
    return dd->size();
}

size_t PricerSolver::get_numberrows_zdd() {
    return zdd->root().row();
}

double PricerSolver::get_cost_edge(int idx) {
    return edges[idx]->cost;
}

void PricerSolver::calculate_edges(scheduleset *set) {
    shared_ptr<node<double>> ptr =
                              zdd_table[zdd->root().row()][zdd->root().col()].list[0];
    size_t count = 0;

    if (set->jobs == NULL) {
        return;
    } else if (!(set->jobs->len)) {
        return;
    }

    Job *tmp_j = (Job *) g_ptr_array_index(set->jobs, count);

    while (ptr->layer != nlayers) {
        job_interval_pair *tmp_pair = (job_interval_pair *) g_ptr_array_index(
                                          ordered_jobs, ptr->layer);
        Job *job = tmp_pair->j;

        if (tmp_j == job) {
            auto e = ptr->out_edge[0].lock();
            g_ptr_array_add(set->e_list, &(e->id));
            count++;

            if (count < set->jobs->len) {
                tmp_j = (Job *) g_ptr_array_index(set->jobs, count);
            }

            ptr = ptr->y;
        } else {
            auto e = ptr->out_edge[1].lock();
            g_ptr_array_add(set->e_list, &(e->id));
            ptr = ptr->n;
        }
    }
}

Optimal_Solution<double> PricerSolver::dynamic_programming_ti(double *pi) {
    Optimal_Solution<double> opt_sol;
    opt_sol.cost = 0;
    double *F;
    Job **A;
    int t_min = 0;
    F = new double [Hmax + 1];
    A = new Job* [Hmax + 1];
    std::vector<Job *> v;


    /** Initialisation */
    F[0] = pi[njobs];
    A[0] = nullptr;

    for (int t = 1; t < Hmax + 1; t++) {
        F[t] = -DBL_MAX / 2;
        A[t] = nullptr;
    }


    /** Recursion */
    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < njobs + 1; i++) {
            int j = i - 1;
            Job *job = (Job *) g_ptr_array_index(jobs, j);

            if (t >=  job->processingime) {
                if (F[t - job->processingime] - (double) value_Fj(t,
                        job) + pi[job->job] >= F[t]) {
                    F[t] = F[t - job->processingime] - value_Fj(t, job) + pi[job->job];
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
        opt_sol.cost += value_Fj(t_min, A[t_min]) ;
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

Optimal_Solution<double> PricerSolver::solve_duration_bdd_double(double *pi) {
    return dd->evaluate_duration(DurationBDDdouble(pi, ordered_jobs, njobs), dd_table);
}

Optimal_Solution<double> PricerSolver::solve_weight_zdd_double(double *pi) {
    return zdd->evaluate_weight(WeightZDDdouble(pi, ordered_jobs, njobs), zdd_table);
}

Optimal_Solution<double> PricerSolver::solve_farkas_double(double *pi) {
    return Optimal_Solution<double>();
}

Optimal_Solution<double> PricerSolver::solve_duration_zdd_double(double *pi){
    return zdd->evaluate_duration(DurationZDDdouble(pi, ordered_jobs, njobs), zdd_duration_table);
}
