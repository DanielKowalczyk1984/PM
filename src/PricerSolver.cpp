#include <PricerSolver.hpp>
#include <PricerConstruct.hpp>
#include <set>
#include <vector>
#include <memory>
/**
 * PricerSolverBase default COnstructor
 */
PricerSolverBase::PricerSolverBase(GPtrArray *_jobs, int _num_machines) :
    jobs(_jobs), njobs(_jobs->len), num_machines(_num_machines), ordered_jobs(nullptr), nlayers(0),
    size_graph(0), nb_removed_edges(0), nb_removed_nodes(0) {
    decision_diagram = nullptr;
}


PricerSolverBase::PricerSolverBase(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    jobs(_jobs), njobs(_jobs->len), num_machines(_num_machines),
    ordered_jobs(_ordered_jobs), nlayers(ordered_jobs->len), size_graph(0),
    nb_removed_edges(0), nb_removed_nodes(0),
    env(new GRBEnv()),
    model(new GRBModel(*env)) {
    /**
     * Construction of decision diagram
     */
    PricerConstruct ps(ordered_jobs);
    decision_diagram = std::unique_ptr<DdStructure<> >(new DdStructure<>(ps));
    size_graph = decision_diagram->size();
}

PricerSolverBase::~PricerSolverBase() {
}

/**
 * Some getters
 */
void PricerSolverBase::iterate_zdd() {
    DdStructure<double>::const_iterator it = decision_diagram->begin();

    for (; it != decision_diagram->end(); ++it) {
        std::set<int>::const_iterator i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << nlayers - *i << " ";
        }

        std::cout << '\n';
    }
}

void PricerSolverBase::print_num_paths() {
    // cout << "Number of paths: " << decision_diagram->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

void PricerSolverBase::create_dot_zdd(const char *name) {
    std::ofstream file;
    file.open(name);
    decision_diagram->dumpDot(file);
    file.close();
}

void PricerSolverBase::print_number_nodes_edges() {
    printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
           nb_removed_nodes);
}

int PricerSolverBase::get_num_remove_nodes() {
    return nb_removed_nodes;
}


int PricerSolverBase::get_num_remove_edges() {
    return nb_removed_edges;
}

size_t PricerSolverBase::get_datasize() {
    return decision_diagram->size();
}

size_t PricerSolverBase::get_size_graph() {
    return size_graph;
}

int PricerSolverBase::get_num_layers() {
    return decision_diagram->topLevel();
}

// size_t PricerSolverBase::get_numberrows_zdd() {
//     return decision_diagram->root().row();
// }

double PricerSolverBase::get_cost_edge(int idx) {
    return 0.0;
}

/**
 * Reduced cost Fixing
 */
void PricerSolverBase::evaluate_nodes(double *pi, int UB, double LB) {
    return;
}

void PricerSolverBase::calculate_edges(scheduleset *set) {
    return;
}

void PricerSolverBase::remove_layers() {
    int first_del = -1;
    int last_del = -1;
    int it = 0;
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
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
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary nodes of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &iter : table[i]) {
            if (!iter.calc_yes) {
                nodeid &cur_node_1 = iter.branch[1];
                cur_node_1 = 0;
            }
        }
    }

    decision_diagram->zddReduce();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram->size();
    printf("The new size of BDD = %lu\n", size_graph);
}

void PricerSolverBase::construct_mipgraph() {
    g.clear();
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
    NodeIdAccessor vertex_nodeid_list(get(boost::vertex_name_t(), g));
    EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), g));

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (size_t j = 0; j < table[i].size(); j++) {
            if (nodeid(i, j) != 0) {
                table[i][j].key = add_vertex(g);
                vertex_nodeid_list[table[i][j].key] = nodeid(i,j);
            }
        }
    }

    int count = 0;
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            if(it.branch[0] != 0) {
                auto &n0 = table.node(it.branch[0]);
                auto a = add_edge(it.key, n0.key,g);
                put(edge_type_list, a.first, false);
                it.low_edge_key = count;
                put(boost::edge_index_t(), g, a.first, count++);
            }

            if(it.branch[1] != 0) {
                auto &n1 = table.node(it.branch[1]);
                auto a = add_edge(it.key, n1.key, g);
                put(edge_type_list, a.first, true);
                it.high_edge_key = count;
                put(boost::edge_index_t(), g, a.first, count++);
            }
        }
    }

    std::cout << "Number of vertices = " << num_vertices(g) << '\n';
    std::cout << "Number of edges = " << num_edges(g) << '\n';
}

void PricerSolverBase::build_mip() {
    try {
        printf("Building Mip model for the extented formulation:\n");
        NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
        IndexAccessor vertex_index_list(get(boost::vertex_index_t(),g));
        NodeIdAccessor vertex_nodeid_list(get(boost::vertex_name_t(),g));
        EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(),g));
        EdgeVarAccessor edge_var_list(get(boost::edge_weight2_t(),g));
        model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
        model->set(GRB_IntParam_Threads, 1);
        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model->set(GRB_IntParam_Presolve, 2);
        model->set(GRB_IntParam_VarBranch, 3);

        /** Constructing variables */
        for (auto it = edges(g); it.first != it.second; it.first++) {
            if (edge_type_list[*it.first]) {
                auto &n = table.node(get(boost::vertex_name_t(),g,source(*it.first, g)));
                double cost = (double) value_Fj(n.GetWeight() + n.GetJob()->processing_time, n.GetJob());
                edge_var_list[*it.first] = model->addVar(0.0, 1.0, cost, GRB_BINARY);
            } else {
                put(boost::edge_weight2_t(),g,*it.first, model->addVar(0.0, (double) num_machines, 0.0, GRB_CONTINUOUS));
            }
        }

        model->update();

        /** Flow constraints */
        size_t num_vertices =  boost::num_vertices(g);
        std::unique_ptr<GRBLinExpr[]> flow_conservation_constr(new GRBLinExpr[num_vertices]());
        std::unique_ptr<char[]> sense_flow(new char[num_vertices]);
        std::unique_ptr<double[]> rhs_flow(new double[num_vertices]);

        for(auto it = vertices(g); it.first != it.second; ++it.first) {
            const auto node_iterator = vertex_nodeid_list[*it.first];
            const auto vertex_index = vertex_index_list[*it.first];
            sense_flow[vertex_index] = GRB_EQUAL;

            auto out_edges_it = boost::out_edges(*it.first, g);
            for(; out_edges_it.first != out_edges_it.second; ++out_edges_it.first) {
                flow_conservation_constr[vertex_index] -= edge_var_list[*out_edges_it.first];
            }

            auto in_edges_it = boost::in_edges(*it.first, g);
            for(; in_edges_it.first != in_edges_it.second; ++in_edges_it.first) {
                flow_conservation_constr[vertex_index] += edge_var_list[*in_edges_it.first];
            }

            if(node_iterator == decision_diagram->root()) {
                rhs_flow[vertex_index] = -(double) num_machines;
            } else if (node_iterator == 1) {
                rhs_flow[vertex_index] = (double) num_machines;
            } else {
                rhs_flow[vertex_index] = 0.0;
            }
        }


        std::unique_ptr<GRBConstr[]> flow_constrs(model->addConstrs(flow_conservation_constr.get(), sense_flow.get(), rhs_flow.get(), nullptr, num_vertices));
        model->update();

        /** Assignment constraints */
        std::unique_ptr<GRBLinExpr[]> assignment(new GRBLinExpr[njobs]());
        std::unique_ptr<char[]> sense(new char[njobs]);
        std::unique_ptr<double[]> rhs(new double[njobs]);

        for (unsigned i = 0; i < jobs->len; ++i) {
            sense[i] = GRB_GREATER_EQUAL;
            rhs[i] = 1.0;
        }

        for (auto it = edges(g); it.first != it.second; it.first++) {
            auto high = edge_type_list[*it.first];
            if (high) {
                auto &n = table.node(get(boost::vertex_name_t(),g,source(*it.first, g)));
                assignment[n.GetJob()->job] += edge_var_list[*it.first];
            }
        }

        std::unique_ptr<GRBConstr[]> assignment_constrs(model->addConstrs(assignment.get(), sense.get(), rhs.get(), nullptr, njobs));
        model->update();

        // model->optimize();
    } catch (GRBException& e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

void PricerSolverBase::calculate_new_ordered_jobs(double *pi, int UB, double LB) {

    /** Remove Layers */
    evaluate_nodes(pi, UB, LB);
    remove_layers();

    /** Construct the new dd (contraction of tables) */
    PricerConstruct ps(ordered_jobs);
    decision_diagram.reset(new DdStructure<>(ps));
    init_table();

    /** Remove nodes for which the high edge points to 0 */
    evaluate_nodes(pi, UB, LB);
    remove_edges();

    construct_mipgraph();
    build_mip();



    // for (int i = njobs  - 1; i >= 0 && count <  8; --i) {
    //     Job *tmp_j = (Job *) g_ptr_array_index(jobs, i);

    //     if (tmp_j->nb_layers == 1) {
    //         zdd->zddSubset(scheduling(tmp_j, ordered_jobs, 2));
    //         zdd->zddReduce();
    //         count++;
    //     }
    // }
}

void PricerSolverBase::construct_lp_sol_from_rmp(const double *columns,
                const GPtrArray *schedule_sets,
                int num_columns,
                double *x) {
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();

    for(int i = 0; i < num_columns; ++i) {
        if(columns[i] > 0.00001) {
            scheduleset *tmp = (scheduleset *) g_ptr_array_index(schedule_sets, i);
            nodeid tmp_nodeid(decision_diagram->root());

            for(unsigned j = 0; j < tmp->job_list->len && tmp_nodeid > 1; ++j) {
                Job* tmp_j = (Job *) g_ptr_array_index(tmp->job_list, j);
                while(true){
                    Node<>& tmp_node = table.node(tmp_nodeid);
                    if(tmp_j == tmp_node.GetJob()) {
                        x[tmp_node.high_edge_key] += columns[i];
                        tmp_nodeid = tmp_node.branch[1];
                        break;
                    } else {
                        x[tmp_node.low_edge_key] += columns[i];
                        tmp_nodeid = tmp_node.branch[0];
                    }
                }
            }
        }
    }
}

bool PricerSolverBase::check_schedule_set(scheduleset *set) {
    int weight = 0;
    int counter = 0;
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
    nodeid tmp_nodeid(decision_diagram->root());

    for(unsigned j = 0; j < set->job_list->len && tmp_nodeid > 1; ++j) {
        Job *tmp_j = (Job *) g_ptr_array_index(set->job_list, j);

        while(true) {
            Node<>& tmp_node = table.node(tmp_nodeid);

            if(tmp_j == tmp_node.GetJob()) {
                tmp_nodeid = tmp_node.branch[1];
                weight += tmp_j->processing_time;
                counter += 1;
                if(weight != table.node(tmp_nodeid).GetWeight() && tmp_nodeid > 1) {
                    return false;
                }
                break;
            } else {
                tmp_nodeid = tmp_node.branch[0];
            }
        }
    }

    return (weight == set->total_processing_time);
}

/**
 * base class for the bdd based pricersolvers
 */
PricerSolverBdd::PricerSolverBdd(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    PricerSolverBase(_jobs, _num_machines, _ordered_jobs) {
    init_table();
    construct_mipgraph();
}

void PricerSolverBdd::init_table() {
    NodeTableEntity<double>& table_new = decision_diagram->getDiagram().privateEntity();

    /** init table */
    table_new.node(decision_diagram->root()).InitNode(0, true);
    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto &it : table_new[i]) {
            if (i != 0) {
                int layer = nlayers - i;
                job_interval_pair *tmp_pair =
                    reinterpret_cast<job_interval_pair *>
                    (g_ptr_array_index(ordered_jobs, layer));
                it.set_job(tmp_pair->j);
                it.set_layer(layer);
                int w = it.GetWeight();
                int p = it.GetJob()->processing_time;

                auto& n0 = table_new.node(it.branch[0]);
                auto& n1 = table_new.node(it.branch[1]);

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
PricerSolverBddSimple::PricerSolverBddSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing BDD with Forward Simple evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = ForwardBddSimpleDouble(njobs);
    reversed_evaluator = BackwardBddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

void PricerSolverBddSimple::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

void PricerSolverBddSimple::evaluate_nodes(double *pi, int UB, double LB) {
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();
            double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];

            if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
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
PricerSolverBddCycle::PricerSolverBddCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing BDD with Forward Cycle evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = ForwardBddCycleDouble(njobs);
    reversed_evaluator = BackwardBddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_forward(evaluator);
}

void PricerSolverBddCycle::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_forward(evaluator);
    decision_diagram->compute_labels_backward(reversed_evaluator);
}

void PricerSolverBddCycle::evaluate_nodes(double *pi, int UB, double LB) {
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();

            if(it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() == job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() == job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
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
PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing BDD with Backward Simple evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = BackwardBddSimpleDouble(njobs);
    reversed_evaluator = ForwardBddSimpleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddBackwardSimple::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverBddBackwardSimple::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardSimple::evaluate_nodes(double *pi, int UB, double LB) {
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();
            double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];

            if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
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
PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
    PricerSolverBdd(_jobs, _num_machines, _ordered_jobs) {
    std::cout << "Constructing BDD with Backward Cycle evaluator" << '\n';
    std::cout << "size BDD = " << get_size_graph() << '\n';
    evaluator = BackwardBddCycleDouble(njobs);
    reversed_evaluator = ForwardBddCycleDouble(njobs);
}

Optimal_Solution<double> PricerSolverBddBackwardCycle::pricing_algorithm(double *_pi) {
    evaluator.initializepi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(double *_pi) {
    evaluator.initializepi(_pi);
    reversed_evaluator.initializepi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

void PricerSolverBddBackwardCycle::evaluate_nodes(double *pi, int UB, double LB) {
    NodeTableEntity<double>& table = decision_diagram->getDiagram().privateEntity();
    compute_labels(pi);
    double reduced_cost = table.node(1).forward_label1.GetF();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto &it : table[i]) {
            int w = it.GetWeight();
            Job *job = it.GetJob();

            if(it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() == job && it.child[1]->backward_label1.get_prev_job() != job) {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label1.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (it.forward_label1.GetPrevJob() != job && it.child[1]->backward_label1.get_prev_job() == job) {
                double result = it.forward_label1.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            } else {
                double result = it.forward_label2.GetF() + it.child[1]->backward_label2.GetF() - value_Fj(w + job->processing_time, job) + pi[job->job] + pi[njobs];
                if (LB - (double)(num_machines - 1)*reduced_cost - result > UB - 1 + 0.0001 && (it.calc_yes)) {
                    it.calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    printf("removed edges = %d\n", nb_removed_edges);
}

/**
 * Pricersolver for the TI index formulation
 */
PricerSolverSimpleDp::PricerSolverSimpleDp(GPtrArray *_jobs, int _num_machines, int _Hmax) :
    PricerSolverBase(_jobs, _num_machines), Hmax(_Hmax), A(new Job* [Hmax + 1]), F(new double[Hmax + 1]) {
        init_table();
}

void PricerSolverSimpleDp::init_table() {
    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < njobs + 1; i++) {
            int j = i - 1;
            Job *job = reinterpret_cast<Job *>(g_ptr_array_index(jobs, j));

            if (t >=  job->processing_time) {
                size_graph++;
            }
        }
    }

    std::cout << "Number of arcs in TI formulation = " << size_graph << '\n';
}

PricerSolverSimpleDp::~PricerSolverSimpleDp() {
}

Optimal_Solution<double> PricerSolverSimpleDp::pricing_algorithm(double *_pi) {
    Optimal_Solution<double> opt_sol;
    opt_sol.cost = 0;
    int t_min = 0;
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

            if (t >=  job->processing_time) {
                if (F[t - job->processing_time]
                    - static_cast<double>(value_Fj(t, job))
                    + _pi[job->job] >= F[t]) {
                    F[t] = F[t - job->processing_time]
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
        t_min -= job->processing_time;
    }

    std::vector<Job *>::reverse_iterator it = v.rbegin();

    for (; it != v.rend(); ++it) {
        g_ptr_array_add(opt_sol.jobs, *it);
    }

    /** Free the memory */
    return opt_sol;
}

/**
 * PricerSolver for the arc-time index formulation
 */
PricerSolverArcTimeDp::PricerSolverArcTimeDp(GPtrArray *_jobs, int _num_machines, int _Hmax) :
    PricerSolverBase(_jobs, _num_machines),
    Hmax(_Hmax),
    n(_jobs->len),
    vector_jobs() {
    for (int i = 0; i < n; ++i) {
        vector_jobs.push_back(reinterpret_cast<Job*>(g_ptr_array_index(jobs, i)));
    }
    job_init(&j0, 0, 0, 0);
    j0.job = n;
    vector_jobs.push_back(&j0);

    init_table();
}


void PricerSolverArcTimeDp::init_table() {
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
        int p = vector_jobs[i]->processing_time;
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
                    && t <= Hmax - tmp->processing_time ) {
                    graph[j][t].insert(it);
                    size_graph++;
                }
            }
        }
    }

    graph[n] = new boost::unordered_set<Job *>[Hmax + 1];
    for (int t = 1; t < Hmax + 1; t++) {
        for (auto &it : vector_jobs) {
            if(t >= it->processing_time) {
                graph[n][t].insert(it);
                size_graph++;
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
            for (int t = tmp_i->processing_time; t <= Hmax - tmp_j->processing_time; ++t) {
                if (delta1(i, j, t) >= 0) {
                    remove_arc(i, j, t);
                    size_graph--;
                } else {
                    remove_arc(j, i, t - tmp_i->processing_time + tmp_j->processing_time);
                    size_graph--;
                }
            }
        }
    }

    for (int j = 0; j < n; ++j) {
        Job *tmp_j = vector_jobs[j];
        for (int t = tmp_j->processing_time; t < Hmax; ++t) {
            if (delta2(j, t) <= 0) {
                remove_arc(n, j, t - tmp_j->processing_time + 1);
                size_graph--;
            } else {
                remove_arc(j, n, t);
                size_graph--;
            }
        }
    }

    for (int j = 0; j < n + 1; ++j) {
        Job* tmp = vector_jobs[j];
        for (int t = 0; t <= Hmax - tmp->processing_time; ++t) {
            if (graph[j][t].empty()) {
                F[j][t] = DBL_MAX/2;
            }
        }
    }
    std::cout << "Number of arcs in ATI formulation = " << size_graph << '\n';
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
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time ) {
                F[j][t] = F[(*it)->job][t - p_matrix[(*it)->job][j]] + value_Fj(t + tmp->processing_time, tmp) - _pi[j];
                A[j][t] = (*it);
                B[j][t] = t - p_matrix[(*it)->job][j];
                it++;
                while(it != graph[j][t].end()) {
                    double result = F[(*it)->job][t - p_matrix[(*it)->job][j]] + value_Fj(t + tmp->processing_time, tmp) - _pi[j];
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
            sol.C_max += vector_jobs[aux_job]->processing_time;
            sol.cost += value_Fj(aux_T + vector_jobs[aux_job]->processing_time, vector_jobs[aux_job]);
            sol.obj += _pi[aux_job] - value_Fj(aux_T + vector_jobs[aux_job]->processing_time, vector_jobs[aux_job]);
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

// /**
//  * Base class for the Zdd based pricing solver
//  */
// PricerSolverZdd::PricerSolverZdd(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
//     PricerSolverBase(_jobs, _num_machines, _ordered_jobs) {
//     init_table();
// }

// ForwardZddNode<double>& PricerSolverZdd::child(tdzdd::NodeId const & id) {
//     return table[id.row()][id.col()];
// }

// void PricerSolverZdd::init_table() {
//     tdzdd::NodeTableHandler<2> &handler = zdd->getDiagram();
//     const tdzdd::NodeTableEntity<2> &diagram = handler.privateEntity();
//     table.init(zdd->topLevel() + 1);

//     /** Init table */
//     for (int i = zdd->topLevel(); i >= 0; i--) {
//         tdzdd::MyVector<tdzdd::Node<2> > const &layer = diagram[i];
//         size_t const m = layer.size();
//         table[i].resize(m);

//         for (auto &it : table[i]) {
//             if (i != 0) {
//                 int layer = nlayers - i;
//                 job_interval_pair *tmp_pair = reinterpret_cast<job_interval_pair *>(
//                     g_ptr_array_index(ordered_jobs, layer));
//                 it.set_job(tmp_pair->j);
//             } else {
//                 it.set_job(nullptr);
//             }
//         }
//     }

//     /**
//      * Init root
//      */
//     table[zdd->topLevel()][0].add_weight(0, 0);

//     /**
//      * Construct in top-down fashion all the nodes in the diagram
//      */
//     for (int i = zdd->topLevel(); i > 0; i--) {
//         size_t const m = table[i].size();

//         for (size_t j = 0; j < m; j++) {
//             Job *job = table[i][j].get_job();
//             tdzdd::NodeId n0 = diagram.child(i, j, 0);
//             tdzdd::NodeId n1 = diagram.child(i, j, 1);

//             for (auto &it : table[i][j].list) {
//                 int w = it->GetWeight();
//                 int p = w + job->processing_time;
//                 it->y = child(n1).add_weight(p, nlayers - n1.row());
//                 it->n = child(n0).add_weight(w, nlayers - n0.row());
//             }
//         }
//     }
// }

// /**
//  * ZDD solver for the flow formulation which ensures that two consecutive jobs are different
//  */
// PricerSolverCycle::PricerSolverCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
//     PricerSolverZdd(_jobs, _num_machines, _ordered_jobs) {
//     evaluator = ForwardZddCycleDouble(njobs);
// }

// Optimal_Solution<double> PricerSolverCycle::pricing_algorithm(double *_pi) {
//     evaluator.initializepi(_pi);
//     return zdd->evaluate_forward(evaluator, table);
// }

// /**
//  * Simple ZDD solver for the flow formulation
//  */
// PricerSolverZddSimple::PricerSolverZddSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs) :
//     PricerSolverZdd(_jobs, _num_machines, _ordered_jobs) {
//     evaluator = ForwardZddSimpleDouble(njobs);
// }

// Optimal_Solution<double> PricerSolverZddSimple::pricing_algorithm(double *_pi) {
//     evaluator.initializepi(_pi);
//     return zdd->evaluate_forward(evaluator, table);
// }

// void PricerSolver::iterate_zdd() {
//     tdzdd::DdStructure<2>::const_iterator it = zdd->begin();

//     for (; it != zdd->end(); ++it) {
//         std::set<int>::const_iterator i = (*it).begin();

//         for (; i != (*it).end(); ++i) {
//             std::cout << nlayers - *i << " ";
//         }

//         std::cout << '\n';
//     }
// }
