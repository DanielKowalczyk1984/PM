#include "PricerSolverBdd.hpp"
#include "PricerConstruct.hpp"
#include "boost/graph/graphviz.hpp"

using namespace std;

PricerSolverBdd::PricerSolverBdd(GPtrArray* _jobs, int _num_machines,
                                 GPtrArray* _ordered_jobs)
    : PricerSolverBase(_jobs, _num_machines, _ordered_jobs),
      size_graph(0),
      nb_removed_edges(0),
      nb_removed_nodes(0),
      env(new GRBEnv()),
      model(new GRBModel(*env))

{
    /**
     * Construction of decision diagram
     */
    PricerConstruct ps(ordered_jobs);
    decision_diagram = std::unique_ptr<DdStructure<>>(new DdStructure<>(ps));
    remove_layers_init();
    decision_diagram->compressBdd();
    size_graph = decision_diagram->size();
    init_table();
    construct_mipgraph();
}

void PricerSolverBdd::construct_mipgraph() {
    mip_graph.clear();
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    NodeIdAccessor   vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (size_t j = 0; j < table[i].size(); j++) {
            if (NodeId(i, j) != 0 && (table[i][j].calc_yes || table[i][j].calc_no)) {
                table[i][j].key = add_vertex(mip_graph);
                vertex_nodeid_list[table[i][j].key] = NodeId(i, j);
            }
        }
    }

    int count = 0;

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (it.branch[0] != 0 && it.calc_no) {
                auto& n0 = table.node(it.branch[0]);
                auto  a = add_edge(it.key, n0.key, mip_graph);
                put(edge_type_list, a.first, false);
                it.low_edge_key = count;
                put(boost::edge_index_t(), mip_graph, a.first, count++);
            }

            if (it.branch[1] != 0 && it.calc_yes) {
                auto& n1 = table.node(it.branch[1]);
                auto  a = add_edge(it.key, n1.key, mip_graph);
                put(edge_type_list, a.first, true);
                it.high_edge_key = count;
                put(boost::edge_index_t(), mip_graph, a.first, count++);
            }
        }
    }

    std::cout << "Number of vertices = " << num_vertices(mip_graph) << '\n';
    std::cout << "Number of edges = " << num_edges(mip_graph) << '\n';
}

void PricerSolverBdd::init_table() {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    /** init table */
    table.node(decision_diagram->root()).init_node(0, true);

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            if (i != 0) {
                int                layer = nlayers - i;
                job_interval_pair* tmp_pair =
                    reinterpret_cast<job_interval_pair*>(
                        g_ptr_array_index(ordered_jobs, layer));
                it.set_job(tmp_pair->j);
                it.set_layer(layer);
                int   w = it.get_weight();
                int   p = it.get_job()->processing_time;
                auto& n0 = table.node(it.branch[0]);
                auto& n1 = table.node(it.branch[1]);
                it.child[0] = n0.init_node(w);
                it.child[1] = n1.init_node(w + p);
            } else {
                it.set_job(nullptr, true);
                it.set_layer(nlayers);
            }
        }
    }
}

void PricerSolverBdd::remove_layers_init() {
    int                first_del = -1;
    int                last_del = -1;
    int                it = 0;
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        if (std::any_of(table[i].begin(), table[i].end(),
                        [](NodeBdd<>& n) { return n.branch[1] != 0; })) {
            if (first_del != -1) {
                g_ptr_array_remove_range(ordered_jobs, first_del,
                                         last_del - first_del + 1);
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
        g_ptr_array_remove_range(ordered_jobs, first_del,
                                 last_del - first_del + 1);
    }

    nlayers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nlayers);
}

void PricerSolverBdd::remove_layers() {
    int                first_del = -1;
    int                last_del = -1;
    int                it = 0;
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        bool remove = true;

        for (auto& iter : table[i]) {
            if (iter.calc_yes) {
                remove = false;
            } else {
                NodeId& cur_node_1 = iter.branch[1];
                cur_node_1 = 0;
            }
        }

        if (!remove) {
            if (first_del != -1) {
                g_ptr_array_remove_range(ordered_jobs, first_del,
                                         last_del - first_del + 1);
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
        g_ptr_array_remove_range(ordered_jobs, first_del,
                                 last_del - first_del + 1);
    }

    nlayers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nlayers);
}

void PricerSolverBdd::remove_edges() {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary nodes of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& iter : table[i]) {
            if (!iter.calc_yes) {
                NodeId& cur_node_1 = iter.branch[1];
                cur_node_1 = 0;
            }

            if (!iter.calc_no)
            {
                NodeId& cur_node_0 = iter.branch[0];
                cur_node_0 = 0;
            }
            
        }
    }

    decision_diagram->compressBdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram->size();
    printf("The new size of BDD = %lu\n", size_graph);
}

void PricerSolverBdd::build_mip() {
    try {
        printf("Building Mip model for the extented formulation:\n");
        NodeTableEntity<>& table =
            decision_diagram->getDiagram().privateEntity();
        IndexAccessor vertex_index_list(
            get(boost::vertex_index_t(), mip_graph));
        NodeIdAccessor vertex_nodeid_list(
            get(boost::vertex_name_t(), mip_graph));
        EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));
        EdgeVarAccessor  edge_var_list(get(boost::edge_weight2_t(), mip_graph));
        model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
        model->set(GRB_IntParam_Threads, 1);
        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model->set(GRB_IntParam_Presolve, 2);
        model->set(GRB_IntParam_VarBranch, 3);

        /** Constructing variables */
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            if (edge_type_list[*it.first]) {
                auto&  n = table.node(get(boost::vertex_name_t(), mip_graph,
                                         source(*it.first, mip_graph)));
                double cost = (double)value_Fj(
                    n.get_weight() + n.get_job()->processing_time, n.get_job());
                edge_var_list[*it.first].x =
                    model->addVar(0.0, 1.0, cost, GRB_CONTINUOUS);
            } else {
                edge_var_list[*it.first].x = model->addVar(
                    0.0, (double)num_machines, 0.0, GRB_CONTINUOUS);
            }
        }

        model->update();
        /** Assignment constraints */
        std::unique_ptr<GRBLinExpr[]> assignment(new GRBLinExpr[njobs]());
        std::unique_ptr<char[]>       sense(new char[njobs]);
        std::unique_ptr<double[]>     rhs(new double[njobs]);

        for (unsigned i = 0; i < jobs->len; ++i) {
            sense[i] = GRB_GREATER_EQUAL;
            rhs[i] = 1.0;
        }

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto high = edge_type_list[*it.first];

            if (high) {
                auto& n = table.node(get(boost::vertex_name_t(), mip_graph,
                                         source(*it.first, mip_graph)));
                assignment[n.get_job()->job] += edge_var_list[*it.first].x;
            }
        }

        std::unique_ptr<GRBConstr[]> assignment_constrs(model->addConstrs(
            assignment.get(), sense.get(), rhs.get(), nullptr, njobs));
        model->update();
        /** Flow constraints */
        size_t num_vertices = boost::num_vertices(mip_graph);
        std::unique_ptr<GRBLinExpr[]> flow_conservation_constr(
            new GRBLinExpr[num_vertices]());
        std::unique_ptr<char[]>   sense_flow(new char[num_vertices]);
        std::unique_ptr<double[]> rhs_flow(new double[num_vertices]);

        for (auto it = vertices(mip_graph); it.first != it.second; ++it.first) {
            const auto node_id = vertex_nodeid_list[*it.first];
            const auto vertex_key = vertex_index_list[*it.first];
            sense_flow[vertex_key] = GRB_EQUAL;
            auto out_edges_it = boost::out_edges(*it.first, mip_graph);

            for (; out_edges_it.first != out_edges_it.second;
                 ++out_edges_it.first) {
                flow_conservation_constr[vertex_key] -=
                    edge_var_list[*out_edges_it.first].x;
            }

            auto in_edges_it = boost::in_edges(*it.first, mip_graph);

            for (; in_edges_it.first != in_edges_it.second;
                 ++in_edges_it.first) {
                flow_conservation_constr[vertex_key] +=
                    edge_var_list[*in_edges_it.first].x;
            }

            if (node_id == decision_diagram->root()) {
                rhs_flow[vertex_key] = -(double)num_machines;
            } else if (node_id == 1) {
                rhs_flow[vertex_key] = (double)num_machines;
            } else {
                rhs_flow[vertex_key] = 0.0;
            }
        }

        std::unique_ptr<GRBConstr[]> flow_constrs(
            model->addConstrs(flow_conservation_constr.get(), sense_flow.get(),
                              rhs_flow.get(), nullptr, num_vertices));
        model->update();
        model->optimize();
        model->write("bdd_formulation.lp");
    } catch (GRBException& e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

void PricerSolverBdd::reduce_cost_fixing(double* pi, int UB, double LB) {
    /** Remove Layers */
    evaluate_nodes(pi, UB, LB);
    remove_layers();
    remove_edges();
    init_table();
    construct_mipgraph();
}

void PricerSolverBdd::add_constraint(Job* job, GPtrArray* list, int order) {
    cout << decision_diagram->size() << '\n';
    scheduling         constr(job, list, order);
    std::ofstream      outf("min1.gv");
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    ColorWriterVertex  vertex_writer(mip_graph, table);
    boost::write_graphviz(outf, mip_graph, vertex_writer);
    decision_diagram->zddSubset(constr);
    outf.close();
    decision_diagram->compressBdd();
    init_table();
    cout << decision_diagram->size() << '\n';
    construct_mipgraph();
    NodeTableEntity<>& table1 = decision_diagram->getDiagram().privateEntity();
    ColorWriterVertex  vertex_writer1(mip_graph, table1);
    outf = std::ofstream("min2.gv");
    boost::write_graphviz(outf, mip_graph, vertex_writer1);
    outf.close();
}

void PricerSolverBdd::construct_lp_sol_from_rmp(const double*    columns,
                                                const GPtrArray* schedule_sets,
                                                int num_columns, double* x) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();

    for (int i = 0; i < num_columns; ++i) {
        if (columns[i] > 0.00001) {
            size_t       counter = 0;
            ScheduleSet* tmp =
                (ScheduleSet*)g_ptr_array_index(schedule_sets, i);
            NodeId tmp_nodeid(decision_diagram->root());

            while (tmp_nodeid > 1) {
                Job* tmp_j;

                if (counter < tmp->job_list->len) {
                    tmp_j = (Job*)g_ptr_array_index(tmp->job_list, counter);
                } else {
                    tmp_j = (Job*)nullptr;
                }

                NodeBdd<>& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    x[tmp_node.high_edge_key] += columns[i];
                    tmp_nodeid = tmp_node.branch[1];
                    counter++;
                } else {
                    x[tmp_node.low_edge_key] += columns[i];
                    tmp_nodeid = tmp_node.branch[0];
                }
            }
        }
    }

    ColorWriterEdge   edge_writer(mip_graph, x);
    ColorWriterVertex vertex_writer(mip_graph, table);
    std::ofstream     outf("min.gv");
    boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    outf.close();
}

double* PricerSolverBdd::project_solution(Solution* sol) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    double*            x = new double[num_edges(mip_graph)]{};

    for (int i = 0; i < sol->nmachines; ++i) {
        size_t     counter = 0;
        GPtrArray* tmp = sol->part[i].machine;
        NodeId     tmp_nodeid(decision_diagram->root());

        while (tmp_nodeid > 1) {
            Job* tmp_j;

            if (counter < tmp->len) {
                tmp_j = (Job*)g_ptr_array_index(tmp, counter);
            } else {
                tmp_j = (Job*)nullptr;
            }

            NodeBdd<>& tmp_node = table.node(tmp_nodeid);

            if (tmp_j == tmp_node.get_job()) {
                x[tmp_node.high_edge_key] += 1.0;
                tmp_nodeid = tmp_node.branch[1];
                counter++;
            } else {
                x[tmp_node.low_edge_key] += 1.0;
                tmp_nodeid = tmp_node.branch[0];
            }
        }
    }

    return x;
}

void PricerSolverBdd::represent_solution(Solution* sol) {
    // double* x = new double[num_edges(g)] {};
    // for(int i = 0; i < sol->nmachines; ++i) {
    //         size_t counter = 0;
    //         GPtrArray *tmp = sol->part[i].machine;
    //         NodeId tmp_nodeid(decision_diagram->root());
    //         while(tmp_nodeid > 1){
    //             Job *tmp_j;
    //             if(counter < tmp->len) {
    //                 tmp_j = (Job *) g_ptr_array_index(tmp, counter);
    //             } else {
    //                 tmp_j = (Job *) nullptr;
    //             }
    //             Node<>& tmp_node = table.node(tmp_nodeid);
    //             if(tmp_j == tmp_node.get_job()) {
    //                 x[tmp_node.high_edge_key] += 1.0;
    //                 tmp_nodeid = tmp_node.branch[1];
    //                 counter++;
    //             } else {
    //                 x[tmp_node.low_edge_key] += 1.0;
    //                 tmp_nodeid = tmp_node.branch[0];
    //             }
    //         }
    // }
    double*            x = project_solution(sol);
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    ColorWriterEdge    edge_writer(mip_graph, x);
    ColorWriterVertex  vertex_writer(mip_graph, table);
    std::ofstream      outf("solution.gv");
    boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    outf.close();
    delete[] x;
}

bool PricerSolverBdd::check_schedule_set(GPtrArray* set) {
    guint              weight = 0;
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    NodeId             tmp_nodeid(decision_diagram->root());

    for (unsigned j = 0; j < set->len && tmp_nodeid > 1; ++j) {
        Job* tmp_j = (Job*)g_ptr_array_index(set, j);

        while (tmp_nodeid > 1) {
            NodeBdd<>& tmp_node = table.node(tmp_nodeid);

            if (tmp_j == tmp_node.get_job()) {
                tmp_nodeid = tmp_node.branch[1];
                weight += 1;

                if (j + 1 != weight) {
                    return false;
                }

                break;
            } else {
                tmp_nodeid = tmp_node.branch[0];
            }
        }
    }

    return (weight == set->len);
}

void PricerSolverBdd::disjunctive_inequality(double* x, Solution* sol) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    int                branch_key = -1;
    std::unique_ptr<GRBModel> model_ineq(new GRBModel(*env));
    EdgeTypeAccessor  edge_type_list(get(boost::edge_weight_t(), mip_graph));
    EdgeVarAccessor   edge_var_list(get(boost::edge_weight2_t(), mip_graph));
    EdgeIndexAccessor edge_index_list(get(boost::edge_index_t(), mip_graph));
    VarsNodeAccessor  node_var_list(get(boost::vertex_distance_t(), mip_graph));
    NodeIdAccessor    node_id_list(get(boost::vertex_name_t(), mip_graph));
    std::unique_ptr<double[]> p(project_solution(sol));
    /**
     * Determine the branch key for the disjunctive program
     */
    int count = 0;

    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto high = edge_type_list[*it.first];
        auto index = edge_index_list[*it.first];

        if (high) {
            if (x[index] > 0.00001 && x[index] < 0.99999 && count < 1) {
                branch_key = index;
                count++;
            }
        }
    }

    printf("branch key = %d\n", branch_key);

    try {
        /**
         * Add variables
         */
        GRBVar s = model_ineq->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        GRBVar t = model_ineq->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            edge_var_list[*it.first].alpha = model_ineq->addVar(
                -GRB_INFINITY, GRB_INFINITY, p[edge_index_list[*it.first]],
                GRB_CONTINUOUS);
        }

        for (auto it = vertices(mip_graph); it.first != it.second; it.first++) {
            node_var_list[*it.first].omega[0] = model_ineq->addVar(
                -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            node_var_list[*it.first].omega[1] = model_ineq->addVar(
                -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        std::unique_ptr<GRBVar[]> pi_0(new GRBVar[njobs]);
        std::unique_ptr<GRBVar[]> pi_1(new GRBVar[njobs]);

        for (int j = 0; j < njobs; j++) {
            pi_0[j] = model_ineq->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0,
                                         GRB_CONTINUOUS);
            pi_1[j] = model_ineq->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0,
                                         GRB_CONTINUOUS);
        }

        GRBVar alpha = model_ineq->addVar(-GRB_INFINITY, GRB_INFINITY, -1.0,
                                          GRB_CONTINUOUS);
        model_ineq->update();
        /**
         * Compute the constraints
         */
        size_t                        num_edges = boost::num_edges(mip_graph);
        std::unique_ptr<GRBLinExpr[]> constraints_0(
            new GRBLinExpr[num_edges + 1]);
        std::unique_ptr<GRBLinExpr[]> constraints_1(
            new GRBLinExpr[num_edges + 1]);
        GRBLinExpr                normalization = GRBLinExpr();
        std::unique_ptr<char[]>   sense(new char[num_edges + 1]);
        std::unique_ptr<double[]> rhs(new double[num_edges + 1]);

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            int  edge_key = edge_index_list[*it.first];
            bool high = edge_type_list[*it.first];
            auto tail = source(*it.first, mip_graph);
            auto head = target(*it.first, mip_graph);
            constraints_0[edge_key] = edge_var_list[*it.first].alpha -
                                      node_var_list[tail].omega[0] +
                                      node_var_list[head].omega[0];
            constraints_1[edge_key] = edge_var_list[*it.first].alpha -
                                      node_var_list[tail].omega[1] +
                                      node_var_list[head].omega[1];

            if (high) {
                NodeId& n = node_id_list[tail];
                constraints_0[edge_key] -= pi_0[table.node(n).get_job()->job];
                constraints_1[edge_key] -= pi_1[table.node(n).get_job()->job];
            }

            if (edge_key == branch_key) {
                constraints_0[edge_key] += s;
                constraints_1[edge_key] -= t;
            }

            normalization += (x[edge_key]) * edge_var_list[*it.first].alpha;
            sense[edge_key] = GRB_GREATER_EQUAL;
            rhs[edge_key] = 0.0;
        }

        normalization -= alpha;
        sense[num_edges] = GRB_GREATER_EQUAL;
        rhs[num_edges] = 0.0;
        int root_key = table.node(decision_diagram->root()).key;
        int terminal_key = table.node(NodeId(1)).key;
        constraints_0[num_edges] += -alpha;
        constraints_1[num_edges] += -alpha;

        for (int j = 0; j < njobs; j++) {
            constraints_0[num_edges] += pi_0[j];
            constraints_1[num_edges] += pi_1[j];
        }

        constraints_0[num_edges] +=
            num_machines * node_var_list[root_key].omega[0] -
            num_machines * node_var_list[terminal_key].omega[0];
        constraints_1[num_edges] +=
            num_machines * node_var_list[root_key].omega[1] -
            num_machines * node_var_list[terminal_key].omega[1] + t;
        /**
         * Add the constraints to the model
         */
        std::unique_ptr<GRBConstr[]> constrs0(
            model_ineq->addConstrs(constraints_0.get(), sense.get(), rhs.get(),
                                   nullptr, num_edges + 1));
        std::unique_ptr<GRBConstr[]> constrs1(
            model_ineq->addConstrs(constraints_1.get(), sense.get(), rhs.get(),
                                   nullptr, num_edges + 1));
        model_ineq->addConstr(normalization, GRB_EQUAL, -1.0);
        model_ineq->update();
        model_ineq->optimize();
        double min = DBL_MAX;

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            double sol = (edge_var_list[*it.first]).alpha.get(GRB_DoubleAttr_X);

            if (CC_ABS(sol) > 0.00001) {
                if (min > CC_ABS(sol)) {
                    min = CC_ABS(sol);
                }
            }
        }

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            double sol = (edge_var_list[*it.first]).alpha.get(GRB_DoubleAttr_X);

            if (CC_ABS(sol) > 0.00001) {
                std::cout << "(" << edge_index_list[*it.first] << ","
                          << sol / min << ") ";
            }
        }

        std::cout << "\n";
        printf("test %f\n", alpha.get(GRB_DoubleAttr_X) / min);
    } catch (GRBException& e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

void PricerSolverBdd::iterate_zdd() {
    DdStructure<NodeBdd<double>>::const_iterator it = decision_diagram->begin();

    for (; it != decision_diagram->end(); ++it) {
        std::set<int>::const_iterator i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << nlayers - *i << " ";
        }

        std::cout << '\n';
    }
}

void PricerSolverBdd::create_dot_zdd(const char* name) {
    std::ofstream file;
    file.open(name);
    decision_diagram->dumpDot(file);
    file.close();
}

void PricerSolverBdd::print_number_nodes_edges() {
    printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
           nb_removed_nodes);
}

int PricerSolverBdd::get_num_remove_nodes() {
    return nb_removed_nodes;
}

int PricerSolverBdd::get_num_remove_edges() {
    return nb_removed_edges;
}

size_t PricerSolverBdd::get_datasize() {
    return decision_diagram->size();
}

size_t PricerSolverBdd::get_size_graph() {
    return decision_diagram->size();
}

int PricerSolverBdd::get_num_layers() {
    return decision_diagram->topLevel();
}

void PricerSolverBdd::print_num_paths() {
    // cout << "Number of paths: " <<
    // decision_diagram->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

double PricerSolverBdd::get_cost_edge(int idx) {
    return 0.0;
}

// void PricerSolverBdd::remove_layers() {
//     int first_del = -1;
//     int last_del = -1;
//     int it = 0;
//     NodeTableEntity<>& table =
//     decision_diagram->getDiagram().privateEntity();

//     /** remove the unnecessary layers of the bdd */
//     for (int i = decision_diagram->topLevel(); i > 0; i--) {
//         bool remove = true;

//         for (auto &iter : table[i]) {
//             if (iter.calc_yes) {
//                 remove = false;
//             } else {
//                 NodeId &cur_node_1 = iter.branch[1];
//                 cur_node_1 = 0;
//             }
//         }

//         if (!remove) {
//             if (first_del != -1) {
//                 g_ptr_array_remove_range(ordered_jobs, first_del, last_del -
//                 first_del + 1); it = it - (last_del - first_del); first_del =
//                 last_del = -1;
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
//         g_ptr_array_remove_range(ordered_jobs, first_del, last_del -
//         first_del + 1);
//     }

//     nlayers = ordered_jobs->len;
//     printf("The new number of layers = %u\n", nlayers);
// }

// void PricerSolverBdd::remove_edges() {
//     NodeTableEntity<>& table =
//     decision_diagram->getDiagram().privateEntity();
//     /** remove the unnecessary nodes of the bdd */
//     for (int i = decision_diagram->topLevel(); i > 0; i--) {
//         for (auto &iter : table[i]) {
//             if (!iter.calc_yes) {
//                 NodeId &cur_node_1 = iter.branch[1];
//                 cur_node_1 = 0;
//             }
//         }
//     }

//     decision_diagram->zddReduce();
//     nb_removed_nodes -= size_graph;
//     size_graph = decision_diagram->size();
//     printf("The new size of BDD = %lu\n", size_graph);
// }

// void PricerSolverBdd::construct_mipgraph() {
//     g.clear();
//     NodeTableEntity<>& table =
//     decision_diagram->getDiagram().privateEntity(); NodeIdAccessor
//     vertex_nodeid_list(get(boost::vertex_name_t(), g)); EdgeTypeAccessor
//     edge_type_list(get(boost::edge_weight_t(), g));

//     for (int i = decision_diagram->topLevel(); i >= 0; i--) {
//         for (size_t j = 0; j < table[i].size(); j++) {
//             if (NodeId(i, j) != 0) {
//                 table[i][j].key = add_vertex(g);
//                 vertex_nodeid_list[table[i][j].key] = NodeId(i, j);
//             }
//         }
//     }

//     int count = 0;
//     for (int i = decision_diagram->topLevel(); i > 0; i--) {
//         for (auto &it : table[i]) {
//             if (it.branch[0] != 0) {
//                 auto &n0 = table.node(it.branch[0]);
//                 auto a = add_edge(it.key, n0.key, g);
//                 put(edge_type_list, a.first, false);
//                 it.low_edge_key = count;
//                 put(boost::edge_index_t(), g, a.first, count++);
//             }

//             if (it.branch[1] != 0) {
//                 auto &n1 = table.node(it.branch[1]);
//                 auto a = add_edge(it.key, n1.key, g);
//                 put(edge_type_list, a.first, true);
//                 it.high_edge_key = count;
//                 put(boost::edge_index_t(), g, a.first, count++);
//             }
//         }
//     }

//     std::cout << "Number of vertices = " << num_vertices(g) << '\n';
//     std::cout << "Number of edges = " << num_edges(g) << '\n';
// }

// void PricerSolverBdd::build_mip() {
//     try {
//         printf("Building Mip model for the extented formulation:\n");
//         NodeTableEntity<>& table =
//         decision_diagram->getDiagram().privateEntity(); IndexAccessor
//         vertex_index_list(get(boost::vertex_index_t(), g)); NodeIdAccessor
//         vertex_nodeid_list(get(boost::vertex_name_t(), g)); EdgeTypeAccessor
//         edge_type_list(get(boost::edge_weight_t(), g)); EdgeVarAccessor
//         edge_var_list(get(boost::edge_weight2_t(), g)); EdgeIndexAccessor
//         edge_index_list(get(boost::edge_index_t(), g));
//         model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
//         model->set(GRB_IntParam_Threads, 1);
//         model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//         model->set(GRB_IntParam_Presolve, 2);
//         model->set(GRB_IntParam_VarBranch, 3);

//         /** Constructing variables */
//         for (auto it = edges(g); it.first != it.second; it.first++) {
//             if (edge_type_list[*it.first]) {
//                 auto &n = table.node(get(boost::vertex_name_t(), g,
//                 source(*it.first, g))); double cost = (double)
//                 value_Fj(n.get_weight() + n.get_job()->processing_time,
//                 n.get_job()); edge_var_list[*it.first].x =
//                 model->addVar(0.0, 1.0, cost, GRB_BINARY);
//             } else {
//                 edge_var_list[*it.first].x = model->addVar(0.0, (double)
//                 num_machines, 0.0, GRB_CONTINUOUS);
//             }
//         }

//         model->update();

//         /** Flow constraints */
//         size_t num_vertices =  boost::num_vertices(g);
//         std::unique_ptr<GRBLinExpr[]> flow_conservation_constr(new
//         GRBLinExpr[num_vertices]()); std::unique_ptr<char[]> sense_flow(new
//         char[num_vertices]); std::unique_ptr<double[]> rhs_flow(new
//         double[num_vertices]);

//         for(auto it = vertices(g); it.first != it.second; ++it.first) {
//             const auto node_id = vertex_nodeid_list[*it.first];
//             const auto vertex_key = vertex_index_list[*it.first];
//             sense_flow[vertex_key] = GRB_EQUAL;

//             auto out_edges_it = boost::out_edges(*it.first, g);
//             for(; out_edges_it.first != out_edges_it.second;
//             ++out_edges_it.first) {
//                 flow_conservation_constr[vertex_key] -=
//                 edge_var_list[*out_edges_it.first].x;
//             }

//             auto in_edges_it = boost::in_edges(*it.first, g);
//             for(; in_edges_it.first != in_edges_it.second;
//             ++in_edges_it.first) {
//                 flow_conservation_constr[vertex_key] +=
//                 edge_var_list[*in_edges_it.first].x;
//             }

//             if(node_id == decision_diagram->root()) {
//                 rhs_flow[vertex_key] = -(double) num_machines;
//             } else if (node_id == 1) {
//                 rhs_flow[vertex_key] = (double) num_machines;
//             } else {
//                 rhs_flow[vertex_key] = 0.0;
//             }
//         }

//         std::unique_ptr<GRBConstr[]>
//         flow_constrs(model->addConstrs(flow_conservation_constr.get(),
//         sense_flow.get(), rhs_flow.get(), nullptr, num_vertices));
//         model->update();

//         /** Assignment constraints */
//         std::unique_ptr<GRBLinExpr[]> assignment(new GRBLinExpr[njobs]());
//         std::unique_ptr<char[]> sense(new char[njobs]);
//         std::unique_ptr<double[]> rhs(new double[njobs]);

//         for (unsigned i = 0; i < jobs->len; ++i) {
//             sense[i] = GRB_GREATER_EQUAL;
//             rhs[i] = 1.0;
//         }

//         for (auto it = edges(g); it.first != it.second; it.first++) {
//             auto high = edge_type_list[*it.first];
//             if (high) {
//                 auto &n =
//                 table.node(get(boost::vertex_name_t(),g,source(*it.first,
//                 g))); assignment[n.get_job()->job] +=
//                 edge_var_list[*it.first].x;
//             }
//         }

//         std::unique_ptr<GRBConstr[]>
//         assignment_constrs(model->addConstrs(assignment.get(), sense.get(),
//         rhs.get(), nullptr, njobs)); model->update();

//         model->optimize();
//         for(auto it = edges(g); it.first != it.second; it.first++) {
//             double sol = (edge_var_list[*it.first]).x.get(GRB_DoubleAttr_X);
//             if(sol > 0.00001) {
//                 printf("test %d %f\n",edge_index_list[*it.first],sol);
//             }
//         }
//     } catch (GRBException& e) {
//         cout << "Error code = " << e.getErrorCode() << endl;
//         cout << e.getMessage() << endl;
//     } catch (...) {
//         cout << "Exception during optimization" << endl;
//     }
// }

// void PricerSolverBdd::reduce_cost_fixing(double *pi, int UB, double LB) {

//     /** Remove Layers */
//     evaluate_nodes(pi, UB, LB);
//     remove_layers();

//     remove_edges();
//     init_table();
//     construct_mipgraph();
// }
