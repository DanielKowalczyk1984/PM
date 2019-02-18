#include <PricerSolver.hpp>
#include <PricerConstruct.hpp>
#include <set>
#include <vector>

/**
 * PricerSolverBase default COnstructor
 */
PricerSolverBase::PricerSolverBase(GPtrArray *_jobs) :
    jobs(_jobs), njobs(_jobs->len), ordered_jobs(nullptr), nlayers(0),
    zdd(nullptr), dd(nullptr), nb_nodes_bdd(0), nb_nodes_zdd(0), nb_arcs_ati(0),
    nb_removed_edges(0), nb_removed_nodes(0) {
    new_dd = nullptr;
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

    new_dd = new DdStructure<>(ps);
    /**
     * Construction of ZDD
     */
    zdd = new tdzdd::DdStructure<2>;
    *zdd = *dd;
    zdd->zddReduce();
    nb_nodes_zdd = zdd->size();
    nb_removed_nodes = 0;
    nb_removed_edges = 0;
    nb_arcs_ati = 0;
}

PricerSolverBase::~PricerSolverBase() {
    if (dd) {
        delete dd;
    }

    if (zdd) {
        delete zdd;
    }

    if(new_dd) {
        delete new_dd;
    }
}

/**
 * Some getters
 */
void PricerSolverBase::IterateZdd() {
    DdStructure<double>::const_iterator it = new_dd->begin();

    for (; it != new_dd->end(); ++it) {
        std::set<int>::const_iterator i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << nlayers - *i << " ";
        }

        std::cout << std::endl;
    }
}

void PricerSolverBase::PrintNumberPaths() {
    cout << "Number of paths: " << dd->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

void PricerSolverBase::create_dot_zdd(const char *name) {
    std::ofstream file;
    file.open(name);
    new_dd->dumpDot(file);
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
    return new_dd->size();
}

size_t PricerSolverBase::get_numberrows_zdd() {
    return new_dd->root().row();
}

int PricerSolverBase::get_nb_arcs_ati() {
    return nb_arcs_ati;
}

double PricerSolverBase::get_cost_edge(int idx) {
    return 0.0;
}

/**
 * Reduced cost Fixing
 */
void PricerSolverBase::evaluate_nodes(double *pi, int UB, double LB, int nmachines) {
    return;
}

// void PricerSolverBase::calculate_new_ordered_jobs() {
//     return;
// }

void PricerSolverBase::calculate_edges(scheduleset *set) {
    return;
}

// void PricerSolverBase::cleanup_table() {
//     NodeTableHandler<double> &handler_bdd = new_dd->getDiagram();
//     NodeTableEntity<double>& table = handler_bdd.privateEntity();
// }

void PricerSolverBase::remove_layers() {
    int first_del = -1;
    int last_del = -1;
    int it = 0;
    NodeTableEntity<double>& table = new_dd->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = new_dd->topLevel(); i > 0; i--) {
        bool remove = true;

        for (auto &iter : table[i]) {
            if (iter.calc_yes) {
                remove = false;
            } else {
                nodeid &cur_node_1 = iter.branch[1];
                cur_node_1 = 0;
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

    nlayers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nlayers);
}

void PricerSolverBase::remove_edges() {
    NodeTableEntity<double>& table_bis = new_dd->getDiagram().privateEntity();

    /** remove the unnecessary nodes of the bdd */
    for (int i = new_dd->topLevel(); i > 0; i--) {
        for (auto &iter : table_bis[i]) {
            if (!iter.calc_yes) {
                nodeid &cur_node_1 = iter.branch[1];
                cur_node_1 = 0;
            }
        }
    }

    new_dd->zddReduce();
    nb_nodes_bdd = new_dd->size();
    printf("The new size of BDD = %lu\n", nb_nodes_bdd);
}

void PricerSolverBase::calculate_new_ordered_jobs(double *pi, int UB, double LB, int nmachines) {
    
    /** Remove Layers */
    evaluate_nodes(pi, UB, LB, nmachines);
    remove_layers();

    /** Construct the new dd (contraction of tables) */
    delete new_dd;
    PricerConstruct ps(ordered_jobs);
    new_dd = new DdStructure<double>(ps);
    InitTable();
    
    /** Remove nodes that for which the high points to 0 */
    evaluate_nodes(pi, UB, LB, nmachines);
    remove_edges();
    


    // for (int i = njobs  - 1; i >= 0 && count <  8; --i) {
    //     Job *tmp_j = (Job *) g_ptr_array_index(jobs, i);

    //     if (tmp_j->nb_layers == 1) {
    //         zdd->zddSubset(scheduling(tmp_j, ordered_jobs, 2));
    //         zdd->zddReduce();
    //         count++;
    //     }
    // }
}

/**
 * base class for the bdd based pricersolvers
 */
PricerSolverBdd::PricerSolverBdd(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBase(_jobs, _ordered_jobs) {
    InitTable();
}

void PricerSolverBdd::InitTable() {
    NodeTableEntity<double>& table_new = new_dd->getDiagram().privateEntity();

    /** init table */
    table_new.node(new_dd->root()).InitNode(0, true);
    for (int i = new_dd->topLevel(); i >= 0; i--) {
        for (auto &it : table_new[i]) {
            if (i != 0) {
                int layer = nlayers - i;
                job_interval_pair *tmp_pair =
                    reinterpret_cast<job_interval_pair *>
                    (g_ptr_array_index(ordered_jobs, layer));
                it.set_job(tmp_pair->j);
                it.set_layer(layer);
                int w = it.GetWeight();
                int p = it.GetJob()->processingime;

                Node<double>& n0 = table_new.node(it.branch[0]);
                Node<double>& n1 = table_new.node(it.branch[1]);

                it.child[0] = n0.InitNode(w);
                it.child[1] = n1.InitNode(w + p);
            } else {
                it.set_job(nullptr, true);
                it.set_layer(nlayers);
            }
        }
    }
}

/**
 *  bdd solver pricersolver for the flow formulation
 */
PricerSolverBddSimple::PricerSolverBddSimple(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _ordered_jobs) {
    std::cout << "Constructing BDD with Forward Simple evaluator" << std::endl;
    std::cout << "size BDD = " << nb_nodes_bdd << " and size ZDD = " << nb_nodes_zdd << std::endl;
    evaluator = ForwardBddSimpleDouble(njobs);
    reversed_evaluator = BackwardBddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return new_dd->evaluate_forward(evaluator);
}

void PricerSolverBddSimple::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    new_dd->compute_labels_forward(evaluator);
    new_dd->compute_labels_backward(reversed_evaluator);
}

void PricerSolverBddSimple::evaluate_nodes(double *pi, int UB, double LB, int nmachines) {
    NodeTableEntity<double>& table = new_dd->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();

    /** check for each node the Lagrangian dual */
    for (int i = new_dd->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();
            double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];

            if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                it.calc_yes = false;
                nb_removed_edges++;
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

/**
 * bdd solver pricersolver for the flow formulation that takes care of the consecutive jobs
 */
PricerSolverBddCycle::PricerSolverBddCycle(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _ordered_jobs) {
    std::cout << "Constructing BDD with Forward Cycle evaluator" << std::endl;
    std::cout << "size BDD = " << nb_nodes_bdd << " and size ZDD = " << nb_nodes_zdd << std::endl;
    evaluator = ForwardBddCycleDouble(njobs);
    reversed_evaluator = BackwardBddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return new_dd->evaluate_forward(evaluator);
}

void PricerSolverBddCycle::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    new_dd->compute_labels_forward(evaluator);
    new_dd->compute_labels_backward(reversed_evaluator);
}

void PricerSolverBddCycle::evaluate_nodes(double *pi, int UB, double LB, int nmachines) {
    NodeTableEntity<double>& table = new_dd->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();

    /** check for each node the Lagrangian dual */
    for (int i = new_dd->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();

            if(it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() == job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() == job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            }

        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

/**
 * backward bdd pricersolver for the flow formulation that takes care of the consecutive jobs
 */
PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _ordered_jobs) {
    std::cout << "Constructing BDD with Backward Simple evaluator" << std::endl;
    std::cout << "size BDD = " << nb_nodes_bdd << " and size ZDD = " << nb_nodes_zdd << std::endl;
    evaluator = BackwardBddSimpleDouble(njobs);
    reversed_evaluator = ForwardBddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddBackwardSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return new_dd->evaluate_backward(evaluator);
}

void PricerSolverBddBackwardSimple::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    new_dd->compute_labels_backward(evaluator);
    new_dd->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double *pi, int UB, double LB, int nmachines) {
    NodeTableEntity<double>& table = new_dd->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();

    /** check for each node the Lagrangian dual */
    for (int i = new_dd->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();
            double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];

            if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                it.calc_yes = false;
                nb_removed_edges++;
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

/**
 * Simple backward bdd pricersolver for the flow formulation
 */
PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _ordered_jobs) {
    std::cout << "Constructing BDD with Backward Cycle evaluator" << std::endl;
    std::cout << "size BDD = " << nb_nodes_bdd << " and size ZDD = " << nb_nodes_zdd << std::endl;
    evaluator = BackwardBddCycleDouble(njobs);
    reversed_evaluator = ForwardBddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddBackwardCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return new_dd->evaluate_backward(evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    new_dd->compute_labels_backward(evaluator);
    new_dd->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double *pi, int UB, double LB, int nmachines) {
    NodeTableEntity<double>& table = new_dd->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();

    /** check for each node the Lagrangian dual */
    for (int i = new_dd->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();

            if(it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() == job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() == job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processingime, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(nmachines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

/**
 * Base class for the Zdd based pricing solver
 */
PricerSolverZdd::PricerSolverZdd(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverBase(_jobs, _ordered_jobs) {
    InitTable();
}

ForwardZddNode<double>& PricerSolverZdd::child(tdzdd::NodeId const & id) {
    return table[id.row()][id.col()];
}

void PricerSolverZdd::InitTable() {
    tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
    const tdzdd::NodeTableEntity<2> &diagram = handler.privateEntity();
    table.init(zdd->topLevel() + 1);

    /** Init table */
    for (int i = zdd->topLevel(); i >= 0; i--) {
        tdzdd::MyVector<tdzdd::Node<2> > const &layer = diagram[i];
        size_t const m = layer.size();
        table[i].resize(m);

        for (auto &it : table[i]) {
            if (i != 0) {
                int layer = nlayers - i;
                job_interval_pair *tmp_pair = reinterpret_cast<job_interval_pair *>(
                    g_ptr_array_index(ordered_jobs, layer));
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
            tdzdd::NodeId n0 = diagram.child(i, j, 0);
            tdzdd::NodeId n1 = diagram.child(i, j, 1);

            for (auto &it : table[i][j].list) {
                int w = it->GetWeight();
                int p = w + job->processingime;
                it->y = child(n1).add_weight(p, nlayers - n1.row());
                it->n = child(n0).add_weight(w, nlayers - n0.row());
            }
        }
    }
}

/**
 * ZDD solver for the flow formulation which ensures that two consecutive jobs are different
 */
PricerSolverCycle::PricerSolverCycle(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverZdd(_jobs, _ordered_jobs) {
    evaluator = ForwardZddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return zdd->evaluate_forward(evaluator, table);
}

/**
 * Simple ZDD solver for the flow formulation
 */
PricerSolverZddSimple::PricerSolverZddSimple(GPtrArray *_jobs, GPtrArray *_ordered_jobs) :
    PricerSolverZdd(_jobs, _ordered_jobs) {
    evaluator = ForwardZddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverZddSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return zdd->evaluate_forward(evaluator, table);
}

/**
 * Pricersolver for the TI index formulation
 */
PricerSolverSimpleDp::PricerSolverSimpleDp(GPtrArray *_jobs, int _Hmax) :
    PricerSolverBase(_jobs), Hmax(_Hmax) {
}

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
    F[0] = -_pi[njobs];
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
                if (F[t - job->processingime]
                    - static_cast<double>(value_Fj(t, job))
                    + _pi[job->job] >= F[t]) {
                    F[t] = F[t - job->processingime]
                           - value_Fj(t, job) + _pi[job->job];
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
        g_ptr_array_add(opt_sol.jobs, *it);
    }

    /** Free the memory */
    delete[] A;
    delete[] F;
    return opt_sol;
}

/**
 * PricerSolver for the arc-time index formulation
 */
PricerSolverArcTimeDp::PricerSolverArcTimeDp(GPtrArray *_jobs, int _Hmax) :
    PricerSolverBase(_jobs),
    Hmax(_Hmax),
    n(_jobs->len),
    vector_jobs() {
    for (int i = 0; i < n; ++i) {
        vector_jobs.push_back(reinterpret_cast<Job*>(g_ptr_array_index(jobs, i)));
    }
    job_init(&j0, 0, 0, 0);
    j0.job = n;
    vector_jobs.push_back(&j0);

    InitTable();
}


void PricerSolverArcTimeDp::InitTable() {
    graph = new boost::unordered_set<Job *>*[n + 1];

    F = new double*[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        F[i] = new double[Hmax + 1]{};
    }

    A = new Job**[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        A[i] = new Job*[Hmax + 1];
    }

    B = new int*[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        B[i] = new int[Hmax + 1];
    }

    p_matrix = new int*[n + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        p_matrix[i] = new int[n + 1];
    }

    for (int i = 0; i < n; ++i) {
        int p = vector_jobs[i]->processingime;
        for (int j = 0; j < n + 1; ++j) {
            p_matrix[i][j] = p;
        }
    }

    for (int j = 0; j < n + 1; ++j) {
        p_matrix[n][j] = (j == n) ? 1 : 0;
    }

    for (int j = 0; j < n; ++j) {
        graph[j] = new boost::unordered_set<Job *>[Hmax + 1];
        Job* tmp = reinterpret_cast<Job*>(g_ptr_array_index(jobs, j));
        for (int t = 0; t < Hmax + 1; t++) {
            for (auto &it : vector_jobs) {
                if (it != tmp
                    && t - p_matrix[it->job][j] >= 0
                    && t <= Hmax - tmp->processingime ) {
                    graph[j][t].insert(it);
                    nb_arcs_ati++;
                }
            }
        }
    }

    graph[n] = new boost::unordered_set<Job *>[Hmax + 1];
    for (int t = 1; t < Hmax + 1; t++) {
        for (auto &it : vector_jobs) {
            if(t >= it->processingime) {
                graph[n][t].insert(it);
                nb_arcs_ati++;
            }
        }
    }

    /**
     * Remove all not needed arcs from the sets
     */
    for (int i = 0; i < n - 1; ++i) {
        Job *tmp_i = vector_jobs[i];
        for (int j = i + 1; j < n; ++j) {
            Job *tmp_j = vector_jobs[j];
            for (int t = tmp_i->processingime; t <= Hmax - tmp_j->processingime; ++t) {
                if (delta1(i, j, t) >= 0) {
                    remove_arc(i, j, t);
                    nb_arcs_ati--;
                } else {
                    remove_arc(j, i, t - tmp_i->processingime + tmp_j->processingime);
                    nb_arcs_ati--;
                }
            }
        }
    }

    for (int j = 0; j < n; ++j) {
        Job *tmp_j = vector_jobs[j];
        for (int t = tmp_j->processingime; t < Hmax; ++t) {
            if (delta2(j, t) <= 0) {
                remove_arc(n, j, t - tmp_j->processingime + 1);
                nb_arcs_ati--;
            } else {
                remove_arc(j, n, t);
                nb_arcs_ati--;
            }
        }
    }

    for (int j = 0; j < n + 1; ++j) {
        Job* tmp = vector_jobs[j];
        for (int t = 0; t <= Hmax - tmp->processingime; ++t) {
            if (graph[j][t].empty()) {
                F[j][t] = DBL_MAX/2;
            }
        }
    }
    std::cout << "Number of arcs in ATI formulation = " << nb_arcs_ati << std::endl;
}

PricerSolverArcTimeDp::~PricerSolverArcTimeDp() {
    for (int i = 0; i < n + 1; ++i) {
        delete[] graph[i];
    }
    delete[] graph;

    for (int i = 0; i < n + 1; ++i) {
        delete[] F[i];
    }
    delete[] F;

    for (int i = 0; i < n + 1; ++i) {
        delete [] A[i];
    }
    delete[] A;

    for (int i = 0; i < n + 1; ++i) {
        delete [] B[i];
    }
    delete[] B;

    for (int i = 0; i < n + 1; ++i) {
        delete [] p_matrix[i];
    }
    delete[] p_matrix;
}

Optimal_Solution<double> PricerSolverArcTimeDp::pricing_algorithm(double *_pi) {
    Optimal_Solution<double> sol(-_pi[n]);
    std::vector<Job*> v;

    F[n][0] = _pi[n];
    double sigma = _pi[n];
    _pi[n] = 0;

    for (int t = 0; t < Hmax + 1; ++t) {
        for (int j = 0; j <= n; ++j) {
            Job *tmp = vector_jobs[j];
            A[j][t] = nullptr;
            B[j][t] = -1;
            F[j][t] = DBL_MAX/2;
            job_iterator it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processingime ) {
                F[j][t] = F[(*it)->job][t - p_matrix[(*it)->job][j]] + value_Fj(t + tmp->processingime, tmp) - _pi[j];
                A[j][t] = (*it);
                B[j][t] = t - p_matrix[(*it)->job][j];
                it++;
                while(it != graph[j][t].end()) {
                    double result = F[(*it)->job][t - p_matrix[(*it)->job][j]] + value_Fj(t + tmp->processingime, tmp) - _pi[j];
                    if (F[j][t] >= result) {
                        F[j][t] = result;
                        A[j][t] = (*it);
                        B[j][t] = t - p_matrix[(*it)->job][j];
                    }
                    it++;
                }
            }
        }
    }

    int job = n;
    int T = Hmax;

    while (T > 0) {
        int aux_job = A[job][T]->job;
        int aux_T = B[job][T];
        if (aux_job != n) {
            v.push_back(vector_jobs[aux_job]);
            sol.C_max += vector_jobs[aux_job]->processingime;
            sol.cost += value_Fj(aux_T + vector_jobs[aux_job]->processingime, vector_jobs[aux_job]);
            sol.obj += _pi[aux_job] - value_Fj(aux_T + vector_jobs[aux_job]->processingime, vector_jobs[aux_job]);
        }
        job = aux_job;
        T = aux_T;
    }


    sol.C_max = 0;
    for (auto &it : v) {
        g_ptr_array_add(sol.jobs, it);
    }
    _pi[n] = sigma;


    return sol;
}

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
