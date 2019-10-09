#include "PricerSolverZdd.hpp"
#include <NodeBddStructure.hpp>
#include "PricerConstruct.hpp"
#include "boost/graph/graphviz.hpp"

PricerSolverZdd::PricerSolverZdd(GPtrArray* _jobs, int _num_machines,
                                 GPtrArray* _ordered_jobs, const char* p_name)
    : PricerSolverBase(_jobs, _num_machines, _ordered_jobs, p_name),
      size_graph(0),
      nb_removed_edges(0),
      nb_removed_nodes(0)

{
    /**
     * Construction of decision diagram
     */
    PricerConstruct ps(ordered_jobs);
    decision_diagram =
        std::unique_ptr<DdStructure<NodeZdd<>>>(new DdStructure<NodeZdd<>>(ps));
    remove_layers_init();
    decision_diagram->compressBdd();
    decision_diagram->reduceZdd();
    size_graph = decision_diagram->size();
    init_table();
    construct_mipgraph();
    lp_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
    solution_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
}

void PricerSolverZdd::construct_mipgraph() {
    mip_graph.clear();
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
    NodeZddIdAccessor vertex_node_zdd_id_list(
        get(boost::vertex_color_t(), mip_graph));
    NodeIdAccessor vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    NodeMipIdAccessor vertex_mip_id_list(
        get(boost::vertex_degree_t(), mip_graph));
    EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (size_t j = 0; j < table[i].size(); j++) {
            auto n{NodeId(i, j)};
            if (n.row() != 0) {
                for (auto& it : table[i][j].list) {
                    auto key = add_vertex(mip_graph);
                    it->key = key;
                    vertex_mip_id_list[it->key] = key;
                    vertex_nodeid_list[it->key] = it->node_id;
                    vertex_node_zdd_id_list[it->key] = it;
                }
            } else {
                if (n != 0) {
                    auto terminal_node = add_vertex(mip_graph);
                    for (auto& it : table[i][j].list) {
                        it->key = terminal_node;
                        vertex_mip_id_list[terminal_node] = terminal_node;
                        vertex_nodeid_list[terminal_node] = it->node_id;
                        vertex_node_zdd_id_list[terminal_node] = it;
                    }
                }
            }
        }
    }

    int count = 0;

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
#ifndef NDEBUG
        auto job = ((job_interval_pair*)g_ptr_array_index(
                        ordered_jobs, ordered_jobs->len - i))
                       ->j;
#endif  // NDEBUG

        for (auto& it : table[i]) {
            if (it.branch[0] != 0) {
                for (auto& iter : it.list) {
                    auto n = iter->n;
                    assert(iter->weight == n->weight);
                    auto a = add_edge(iter->key, n->key, mip_graph);
                    put(edge_type_list, a.first, false);
                    iter->low_edge_key = count;
                    put(boost::edge_index_t(), mip_graph, a.first, count++);
                }
            }

            if (it.branch[1] != 0) {
                for (auto& iter : it.list) {
                    auto y = iter->y;
                    assert(iter->weight + job->processing_time == y->weight);
                    auto a = add_edge(iter->key, y->key, mip_graph);
                    put(edge_type_list, a.first, true);
                    iter->high_edge_key = count;
                    put(boost::edge_index_t(), mip_graph, a.first, count++);
                }
            }
        }
    }

    std::cout << "Number of vertices = " << num_vertices(mip_graph) << '\n';
    std::cout << "Number of edges = " << num_edges(mip_graph) << '\n';
}

void PricerSolverZdd::init_table() {
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
    /** init table */
    auto& n = table.node(decision_diagram->root());
    n.add_sub_node(0, decision_diagram->root(), true, false);
    n.set_node_id(decision_diagram->root());
    if (table.node(decision_diagram->root()).list.size() > 1) {
        table.node(decision_diagram->root()).list.pop_back();
    }

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        int                layer = nb_layers - i;
        job_interval_pair* tmp_pair = reinterpret_cast<job_interval_pair*>(
            g_ptr_array_index(ordered_jobs, layer));

        for (auto& it : table[i]) {
            if (i != 0) {
                it.set_job(tmp_pair->j);
                it.set_layer(layer);
                auto& n0 = table.node(it.branch[0]);
                auto& n1 = table.node(it.branch[1]);
                int   p = it.get_job()->processing_time;
                it.child[0] = table.node_ptr(it.branch[0]);
                it.child[1] = table.node_ptr(it.branch[1]);
                for (auto& iter : it.list) {
                    int w = iter->weight;
                    iter->n = n0.add_weight(w, it.branch[0]);
                    iter->y = n1.add_weight(w + p, it.branch[1]);
                }
            } else {
                it.set_job(nullptr, true);
                it.set_layer(nb_layers);
                it.set_root_node(true);
            }
        }
    }
}

void PricerSolverZdd::remove_layers_init() {
    int                         first_del = -1;
    int                         last_del = -1;
    int                         it = 0;
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        if (std::any_of(table[i].begin(), table[i].end(),
                        [](NodeZdd<>& n) { return n.branch[1] != 0; })) {
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

    nb_layers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nb_layers);
}

void PricerSolverZdd::remove_layers() {
    int                         first_del = -1;
    int                         last_del = -1;
    int                         it = 0;
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        bool remove_layer = true;

        for (auto& iter : table[i]) {
            auto end = std::remove_if(
                iter.list.begin(), iter.list.end(),
                [](std::shared_ptr<SubNodeZdd<>> n) { return !(n->calc_yes); });
            iter.list.erase(end, iter.list.end());

            if (iter.list.empty()) {
                NodeId& cur_node_1 = iter.branch[1];
                cur_node_1 = 0;
            } else {
                remove_layer = false;
            }
        }

        if (!remove_layer) {
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

    nb_layers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nb_layers);
}

void PricerSolverZdd::remove_edges() {
    decision_diagram->compressBdd();
    decision_diagram->reduceZdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram->size();
    printf("The new size of ZDD = %lu\n", size_graph);
}

void PricerSolverZdd::build_mip() {
    try {
        printf("Building Mip model for the extended formulation:\n");
        NodeIdAccessor vertex_nodeid_list(
            get(boost::vertex_name_t(), mip_graph));
        NodeMipIdAccessor vertex_mip_id_list(
            get(boost::vertex_degree_t(), mip_graph));
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
                auto& n = get(boost::vertex_color_t(), mip_graph,
                              source(*it.first, mip_graph));
                Job*  job = n->get_job();

                double cost =
                    (double)value_Fj(n->weight + job->processing_time, job);
                edge_var_list[*it.first].x =
                    model->addVar(0.0, 1.0, cost, GRB_CONTINUOUS);
            } else {
                edge_var_list[*it.first].x = model->addVar(
                    0.0, (double)num_machines, 0.0, GRB_CONTINUOUS);
            }
        }

        model->update();
        /** Assignment constraints */
        std::unique_ptr<GRBLinExpr[]> assignment(new GRBLinExpr[nb_jobs]());
        std::unique_ptr<char[]>       sense(new char[nb_jobs]);
        std::unique_ptr<double[]>     rhs(new double[nb_jobs]);

        for (unsigned i = 0; i < jobs->len; ++i) {
            sense[i] = GRB_GREATER_EQUAL;
            rhs[i] = 1.0;
        }

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto high = edge_type_list[*it.first];

            if (high) {
                auto& n = get(boost::vertex_color_t(), mip_graph,
                              source(*it.first, mip_graph));
                assignment[n->get_job()->job] += edge_var_list[*it.first].x;
            }
        }

        std::unique_ptr<GRBConstr[]> assignment_constrs(model->addConstrs(
            assignment.get(), sense.get(), rhs.get(), nullptr, nb_jobs));
        model->update();
        /** Flow constraints */
        size_t num_vertices = boost::num_vertices(mip_graph);
        // boost::num_vertices(mip_graph) - table[0][1].list.size() + 1;
        std::unique_ptr<GRBLinExpr[]> flow_conservation_constr(
            new GRBLinExpr[num_vertices]());
        std::unique_ptr<char[]>   sense_flow(new char[num_vertices]);
        std::unique_ptr<double[]> rhs_flow(new double[num_vertices]);

        for (auto it = vertices(mip_graph); it.first != it.second; ++it.first) {
            const auto node_id = vertex_nodeid_list[*it.first];
            const auto vertex_key = vertex_mip_id_list[*it.first];
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
            } else if (node_id.row() == 0) {
                rhs_flow[vertex_key] = (double)num_machines;
            } else {
                rhs_flow[vertex_key] = 0.0;
            }
        }

        std::unique_ptr<GRBConstr[]> flow_constrs(
            model->addConstrs(flow_conservation_constr.get(), sense_flow.get(),
                              rhs_flow.get(), nullptr, num_vertices));
        model->update();
        model->write("zdd_formulation.lp");
        model->optimize();
    } catch (GRBException& e) {
        std::cout << "Error code = " << e.getErrorCode() << "\n";
        std::cout << e.getMessage() << "\n";
    } catch (...) {
        std::cout << "Exception during optimization"
                  << "\n";
    }
}

void PricerSolverZdd::reduce_cost_fixing(double* pi, int UB, double LB) {
    /** Remove Layers */
    evaluate_nodes(pi, UB, LB);
    remove_layers();
    remove_edges();
    init_table();
    evaluate_nodes(pi, UB, LB);
    remove_edges();
    init_table();
    construct_mipgraph();
}

void PricerSolverZdd::add_constraint(Job* job, GPtrArray* list, int order) {
    std::cout << decision_diagram->size() << '\n';
    scheduling constr(job, list, order);
    // std::ofstream outf("min1.gv");
    // NodeTableEntity<NodeZdd<>>& table =
    // decision_diagram->getDiagram().privateEntity(); ColorWriterVertex
    // vertex_writer(g, table); boost::write_graphviz(outf, g, vertex_writer);
    decision_diagram->zddSubset(constr);
    // outf.close();
    decision_diagram->compressBdd();
    init_table();
    std::cout << decision_diagram->size() << '\n';
    construct_mipgraph();
    // NodeTableEntity<NodeZdd<>>& table1 =
    // decision_diagram->getDiagram().privateEntity(); ColorWriterVertex
    // vertex_writer1(g, table1); outf = std::ofstream("min2.gv");
    // boost::write_graphviz(outf, g, vertex_writer1);
    // outf.close();
}

void PricerSolverZdd::construct_lp_sol_from_rmp(const double*    columns,
                                                const GPtrArray* schedule_sets,
                                                int              num_columns) {
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
    std::fill(lp_x.get(), lp_x.get() + get_nb_edges(), 0.0);
    for (int i = 0; i < num_columns; ++i) {
        if (columns[i] > 0.00001) {
            size_t       counter = 0;
            ScheduleSet* tmp =
                (ScheduleSet*)g_ptr_array_index(schedule_sets, i);
            NodeId                        tmp_nodeid(decision_diagram->root());
            std::shared_ptr<SubNodeZdd<>> tmp_sub_node =
                table.node(tmp_nodeid).list[0];

            while (tmp_nodeid > 1) {
                Job* tmp_j;

                if (counter < tmp->job_list->len) {
                    tmp_j = (Job*)g_ptr_array_index(tmp->job_list, counter);
                } else {
                    tmp_j = (Job*)nullptr;
                }

                NodeZdd<>& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    lp_x[tmp_sub_node->high_edge_key] += columns[i];
                    tmp_nodeid = tmp_node.branch[1];
                    tmp_sub_node = tmp_sub_node->y;
                    counter++;
                } else {
                    lp_x[tmp_sub_node->low_edge_key] += columns[i];
                    tmp_nodeid = tmp_node.branch[0];
                    tmp_sub_node = tmp_sub_node->n;
                }
            }
        }
    }

    // ColorWriterEdge edge_writer(g, x);
    // ColorWriterVertex vertex_writer(g, table);
    // std::ofstream outf("min.gv");
    // boost::write_graphviz(outf, g, vertex_writer, edge_writer);
    // outf.close();
}

void PricerSolverZdd::project_solution(Solution* sol) {
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
    // double* x = new double[num_edges(mip_graph)]{};
    std::fill(solution_x.get(), solution_x.get() + get_nb_edges(), 0.0);

    for (int i = 0; i < sol->nb_machines; ++i) {
        size_t                        counter = 0;
        GPtrArray*                    tmp = sol->part[i].machine;
        NodeId                        tmp_nodeid(decision_diagram->root());
        std::shared_ptr<SubNodeZdd<>> tmp_sub_node =
            table.node(tmp_nodeid).list[0];

        while (tmp_nodeid > 1) {
            Job* tmp_j;

            if (counter < tmp->len) {
                tmp_j = (Job*)g_ptr_array_index(tmp, counter);
            } else {
                tmp_j = (Job*)nullptr;
            }

            NodeZdd<>& tmp_node = table.node(tmp_nodeid);

            if (tmp_j == tmp_node.get_job()) {
                solution_x[tmp_sub_node->high_edge_key] += 1.0;
                tmp_nodeid = tmp_node.branch[1];
                tmp_sub_node = tmp_sub_node->y;
                counter++;
            } else {
                solution_x[tmp_sub_node->low_edge_key] += 1.0;
                tmp_nodeid = tmp_node.branch[0];
                tmp_sub_node = tmp_sub_node->n;
            }
        }
    }
}

void PricerSolverZdd::represent_solution(Solution* sol) {
    // double* x = new double[num_edges(g)] {};
    // for(int i = 0; i < sol->nb_machines; ++i) {
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
    project_solution(sol);
    // NodeTableEntity<NodeZdd<>>& table =
    // decision_diagram->getDiagram().privateEntity(); ColorWriterEdge
    // edge_writer(g, x); ColorWriterVertex vertex_writer(g, table);
    // std::ofstream outf("solution.gv"); boost::write_graphviz(outf, g,
    // vertex_writer, edge_writer); outf.close();
    // delete[] x;
}

bool PricerSolverZdd::check_schedule_set(GPtrArray* set) {
    guint                       weight = 0;
    NodeTableEntity<NodeZdd<>>& table =
        decision_diagram->getDiagram().privateEntity();
    NodeId tmp_nodeid(decision_diagram->root());

    for (unsigned j = 0; j < set->len; ++j) {
        Job* tmp_j = (Job*)g_ptr_array_index(set, j);
        while (tmp_nodeid > 1) {
            NodeZdd<>& tmp_node = table.node(tmp_nodeid);

            if (tmp_j == tmp_node.get_job()) {
                tmp_nodeid = tmp_node.branch[1];
                weight += 1;

                if (j + 1 != weight) {
                    return false;
                }
            } else {
                tmp_nodeid = tmp_node.branch[0];
            }
        }
    }

    return (weight == set->len);
}

void PricerSolverZdd::disjunctive_inequality(double* x, Solution* sol) {}

void PricerSolverZdd::iterate_zdd() {
    DdStructure<NodeZdd<double>>::const_iterator it = decision_diagram->begin();

    for (; it != decision_diagram->end(); ++it) {
        std::set<int>::const_iterator i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << nb_layers - *i << " ";
        }

        std::cout << '\n';
    }
}

void PricerSolverZdd::create_dot_zdd(const char* name) {
    std::ofstream file;
    file.open(name);
    decision_diagram->dumpDot(file);
    file.close();
}

void PricerSolverZdd::print_number_nodes_edges() {
    printf("removed edges = %d, removed nodes = %d\n", nb_removed_edges,
           nb_removed_nodes);
}

int PricerSolverZdd::get_num_remove_nodes() {
    return nb_removed_nodes;
}

int PricerSolverZdd::get_num_remove_edges() {
    return nb_removed_edges;
}

size_t PricerSolverZdd::get_nb_edges() {
    return num_edges(mip_graph);
}

size_t PricerSolverZdd::get_nb_vertices() {
    return decision_diagram->size();
}

int PricerSolverZdd::get_num_layers() {
    return decision_diagram->topLevel();
}

void PricerSolverZdd::print_num_paths() {}

int PricerSolverZdd::get_int_attr_model(enum MIP_Attr c) {
    int val = -1;
    switch (c) {
        case MIP_Attr_Nb_Vars:
            val = model->get(GRB_IntAttr_NumVars);
            break;
        case MIP_Attr_Nb_Constr:
            val = model->get(GRB_IntAttr_NumConstrs);
            break;
        case MIP_Attr_Status:
            val = model->get(GRB_IntAttr_Status);
            break;
        default:
            break;
    }

    return val;
}

double PricerSolverZdd::get_dbl_attr_model(enum MIP_Attr c) {
    double val = -1.0;
    switch (c) {
        case MIP_Attr_Obj_Bound:
            val = model->get(GRB_DoubleAttr_ObjBound);
            break;
        case MIP_Attr_Obj_Bound_LP:
            val = model->get(GRB_DoubleAttr_ObjBoundC);
            break;
        case MIP_Attr_Mip_Gap:
            val = model->get(GRB_DoubleAttr_MIPGap);
            break;
        case MIP_Attr_Run_Time:
            val = model->get(GRB_DoubleAttr_Runtime);
            break;
        case MIP_Attr_Nb_Simplex_Iter:
            val = model->get(GRB_DoubleAttr_IterCount);
            break;
        case MIP_Attr_Nb_Nodes:
            val = model->get(GRB_DoubleAttr_NodeCount);
            break;
        default:
            break;
    }

    return val;
}
