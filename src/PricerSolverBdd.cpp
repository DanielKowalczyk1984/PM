#include "PricerSolverBdd.hpp"
#include <list>
#include "PricerConstruct.hpp"
#include "boost/graph/graphviz.hpp"

using namespace std;

PricerSolverBdd::PricerSolverBdd(GPtrArray* _jobs, int _num_machines,
                                 GPtrArray* _ordered_jobs, const char* p_name)
    : PricerSolverBase(_jobs, _num_machines, _ordered_jobs, p_name),
      size_graph(0),
      nb_removed_edges(0),
      nb_removed_nodes(0),
      H_min(0)

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
    calculate_H_min();
    cleanup_arcs();
    topdown_filtering();
    check_infeasible_arcs();
    bottum_up_filtering();
    construct_mipgraph();
    lp_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
    solution_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
}

PricerSolverBdd::PricerSolverBdd(GPtrArray* _jobs, int _nb_machines,
                                 GPtrArray* _ordered_jobs, int* _take_jobs,
                                 int _Hmax, const char* p_name)
    : PricerSolverBase(_jobs, _nb_machines, _ordered_jobs, p_name) {
    PricerConstructTI ps(ordered_jobs, _take_jobs, _Hmax);
    decision_diagram = std::unique_ptr<DdStructure<>>(new DdStructure<>(ps));
    remove_layers_init();
    decision_diagram->compressBdd();
    init_table();
    calculate_H_min();
    cleanup_arcs();
    topdown_filtering();
    check_infeasible_arcs();
    bottum_up_filtering();
    construct_mipgraph();
    lp_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
    solution_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
    H_min = 0;
}
int g_compare_duration(gconstpointer a, gconstpointer b) {
    const Job* x = *((Job* const*)a);
    const Job* y = *((Job* const*)b);

    if (x->processing_time < y->processing_time) {
        return -1;
    } else {
        return 1;
    }
}

void PricerSolverBdd::calculate_H_min() {
    int        p_sum = 0;
    GPtrArray* duration = g_ptr_array_new();
    for (int j = 0; j < nb_jobs; j++) {
        Job* job = (Job*)g_ptr_array_index(jobs, j);
        g_ptr_array_add(duration, job);
        p_sum += job->processing_time;
    }
    g_ptr_array_sort(duration, g_compare_duration);

    int    m = 0;
    int    i = nb_jobs - 1;
    double tmp = p_sum;
    do {
        Job* job = (Job*)g_ptr_array_index(duration, nb_jobs - 1);
        tmp -= job->processing_time;
        std::cout << job->processing_time << " " << m << "\n";
        m++;
    } while (m < num_machines - 1);

    H_min = (int)floor(tmp / num_machines);

    g_ptr_array_free(duration, TRUE);
}

void PricerSolverBdd::construct_mipgraph() {
    mip_graph.clear();
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    NodeIdAccessor   vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (size_t j = 0; j < table[i].size(); j++) {
            if (NodeId(i, j) != 0
                // && (table[i][j].calc_yes || table[i][j].calc_no)
            ) {
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
    NodeBdd<>& root = table.node(decision_diagram->root());
    root.init_node(0, true);
    root.all = boost::dynamic_bitset<>{nb_jobs, 0};

    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            if (i != 0) {
                int                layer = nb_layers - i;
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
                it.set_layer(nb_layers);
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

    nb_layers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nb_layers);
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

    nb_layers = ordered_jobs->len;
    printf("The new number of layers = %u\n", nb_layers);
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

            if (!iter.calc_no) {
                NodeId& cur_node_0 = iter.branch[0];
                cur_node_0 = 0;
            }
        }
    }

    decision_diagram->compressBdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram->size();
    printf("The new size of BDD = %lu\n", size_graph);
    std::cout
        << "-------------------------------------------------------------\n";
}

void PricerSolverBdd::print_representation_file() {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    NodeIdAccessor   vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));
    EdgeIndexAccessor edge_index_list(get(boost::edge_index_t(), mip_graph));
    std::unique_ptr<vector<int>[]> index_edge(new vector<int>[nb_jobs]);

    string outfile_file_mip_str =
        problem_name + "_" + std::to_string(num_machines) + ".txt";
    std::ofstream out_file_mip(outfile_file_mip_str);

    out_file_mip << boost::num_vertices(mip_graph) << " "
                 << boost::num_edges(mip_graph) << " " << nb_jobs << " "
                 << num_machines << "\n\n";

    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto& head = table.node(get(boost::vertex_name_t(), mip_graph,
                                    source(*it.first, mip_graph)));
        auto& n = table.node(get(boost::vertex_name_t(), mip_graph,
                                 target(*it.first, mip_graph)));
        auto  high = edge_type_list[*it.first];
        if (high) {
            auto cost = (double)value_Fj(
                head.get_weight() + head.get_job()->processing_time,
                head.get_job());
            out_file_mip << head.key << " " << n.key << " " << cost << "\n";
            index_edge[n.get_job()->job].push_back(edge_index_list[*it.first]);
        } else {
            out_file_mip << head.key << " " << n.key << " " << 0.0 << "\n";
        }
    }

    out_file_mip << "\n";

    for (int i = 0; i < nb_jobs; i++) {
        for (auto& it : index_edge[i]) {
            out_file_mip << it << " ";
        }
        out_file_mip << "\n";
    }

    out_file_mip << "\n";

    for (auto it = vertices(mip_graph); it.first != it.second; it.first++) {
        const auto node_id = vertex_nodeid_list[*it.first];
        auto&      n = table.node(node_id);
        if (node_id > 1) {
            out_file_mip << n.get_job()->job << " " << n.get_weight() << "\n";
        }
    }

    out_file_mip.close();
}

void PricerSolverBdd::build_mip() {
    try {
        printf("Building Mip model for the extended formulation:\n");
        NodeTableEntity<>& table =
            decision_diagram->getDiagram().privateEntity();
        IndexAccessor vertex_index_list(
            get(boost::vertex_index_t(), mip_graph));
        NodeIdAccessor vertex_nodeid_list(
            get(boost::vertex_name_t(), mip_graph));
        EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));
        EdgeVarAccessor  edge_var_list(get(boost::edge_weight2_t(), mip_graph));
        EdgeIndexAccessor edge_index_list(
            get(boost::edge_index_t(), mip_graph));

        /** Constructing variables */
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            if (edge_type_list[*it.first]) {
                auto&  n = table.node(get(boost::vertex_name_t(), mip_graph,
                                         source(*it.first, mip_graph)));
                auto   C = n.get_weight() + n.get_job()->processing_time;
                double cost = (double)value_Fj(C, n.get_job());
                edge_var_list[*it.first].x =
                    model->addVar(0.0, 1.0, cost, GRB_BINARY);
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
            sense[i] = GRB_EQUAL;
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
            assignment.get(), sense.get(), rhs.get(), nullptr, nb_jobs));
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
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            // edge_var_list[*it.first].x.set(GRB_DoubleAttr_PStart,
            //                                lp_x[edge_index_list[*it.first]]);
            edge_var_list[*it.first].x.set(
                GRB_DoubleAttr_Start, solution_x[edge_index_list[*it.first]]);
        }
        model->write(problem_name + "_" + std::to_string(num_machines) +
                     "_mip.mps");
        model->optimize();

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            int index = edge_index_list[*it.first];
            solution_x[index] =
                edge_var_list[*it.first].x.get(GRB_DoubleAttr_X);
        }

        ColorWriterEdgeX  edge_writer(mip_graph, solution_x.get());
        ColorWriterVertex vertex_writer(mip_graph, table);
        string            file_name = "lp_solution_" + problem_name + "_" +
                           std::to_string(num_machines) + ".gv";
        std::ofstream outf(file_name);
        boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
        outf.close();

        ColorWriterEdgeIndex edge_writer_index(mip_graph);
        file_name = "index_" + problem_name + "_" +
                    std::to_string(num_machines) + ".gv";
        std::ofstream outf_index(file_name);
        boost::write_graphviz(outf_index, mip_graph, vertex_writer,
                              edge_writer_index);
        outf_index.close();
    } catch (GRBException& e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}

void PricerSolverBdd::reduce_cost_fixing(double* pi, int UB, double LB) {
    /** Remove Layers */
    std::cout << "Starting Reduced cost fixing\n";
    evaluate_nodes(pi, UB, LB);
    topdown_filtering();
    check_infeasible_arcs();
    // bottumup_filtering();
    cleanup_arcs();

    construct_mipgraph();
    // equivalent_paths_filtering();
}

void PricerSolverBdd::cleanup_arcs() {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();

    table.node(0).backward_distance[0] = INT_MIN;
    table.node(0).backward_distance[1] = INT_MIN;
    table.node(1).backward_distance[0] = 0;
    table.node(1).backward_distance[1] = 0;
    bool removed_edges = false;
    int  nb_edges_removed_tmp = 0;

    for (int i = 1; i <= decision_diagram->topLevel(); i++) {
        for (auto& it : table[i]) {
            it.calc_no = true;
            it.calc_yes = true;
            NodeBdd<>& cur_node_0 = table.node(it.branch[0]);
            NodeBdd<>& cur_node_1 = table.node(it.branch[1]);

            if (cur_node_0.backward_distance[0] <
                cur_node_0.backward_distance[1]) {
                it.backward_distance[0] = cur_node_0.backward_distance[1];
            } else {
                it.backward_distance[0] = cur_node_0.backward_distance[0];
            }

            int result0 =
                cur_node_1.backward_distance[0] + it.get_job()->processing_time;
            int result1 =
                cur_node_1.backward_distance[1] + it.get_job()->processing_time;

            if (result0 < result1) {
                it.backward_distance[1] = result1;
            } else {
                it.backward_distance[1] = result0;
            }
        }
    }
    /** remove the unnecessary nodes of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& iter : table[i]) {
            if (iter.get_weight() + iter.backward_distance[0] < H_min) {
                iter.calc_no = false;
                removed_edges = true;
                nb_edges_removed_tmp++;
                nb_removed_edges++;
                // iter.calc_yes = false;
                // removed_edges = true;
                // nb_edges_removed_tmp++;
                // nb_removed_edges++;
            }

            // if (iter.get_weight() + iter.backward_distance[1] < H_min &&
            // iter.branch[1] == 1) {
            //     iter.calc_yes = false;
            //     removed_edges = true;
            //     nb_edges_removed_tmp++;
            //     nb_removed_edges++;
            // }
        }
    }

    if (removed_edges) {
        std::cout << "Number of edges removed by cleanup arcs = "
                  << nb_edges_removed_tmp << "\n";
        std::cout << "Number of edges removed in total = " << nb_removed_edges
                  << "\n";
        remove_layers();
        remove_edges();
        init_table();
    }
}

void PricerSolverBdd::topdown_filtering() {
    bool               removed_edges = false;
    int                nb_edges_removed_tmp = 0;
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    NodeBdd<>&         root = table.node(decision_diagram->root());
    root.init_node(0, true);
    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{nb_jobs, 0};
            it.calc_yes = true;
        }
    }

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto& n0 = table.node(it.branch[0]);
            if (n0.visited) {
                n0.all &= it.all;
            } else {
                n0.visited = true;
                n0.all = boost::dynamic_bitset<>{nb_jobs, 0};
                n0.all |= it.all;
            }
            auto& n1 = table.node(it.branch[1]);
            if (n1.visited) {
                if (n1.all[it.get_job()->job]) {
                    n1.all &= it.all;
                    n1.all[it.get_job()->job] = 1;
                } else {
                    n1.all &= it.all;
                }
            } else {
                n1.all = boost::dynamic_bitset<>{nb_jobs, 0};
                n1.all |= it.all;
                n1.all[it.get_job()->job] = 1;
                n1.visited = true;
            }
        }
    }

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (it.all[it.get_job()->job]) {
                removed_edges = true;
                it.calc_yes = false;
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        std::cout << "removing edges based on top-down iteration\n";
        std::cout << "Number edges removed top-bottom = "
                  << nb_edges_removed_tmp << "\n";
        std::cout << "Number edges removed total = " << nb_removed_edges
                  << "\n";
        remove_layers();
        remove_edges();
        init_table();
        cleanup_arcs();
        // continue;
    }
}

void PricerSolverBdd::bottum_up_filtering() {
    bool              removed_edges = false;
    int               nb_edges_removed_tmp = 0;
    NodeTableEntity<> table = decision_diagram->getDiagram().privateEntity();
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{nb_jobs, 0};
            it.calc_yes = true;
        }
    }

    table.node(0).all = boost::dynamic_bitset<>{nb_jobs, 0};
    table.node(1).all = boost::dynamic_bitset<>{nb_jobs, 0};
    table.node(0).all.flip();

    for (int i = 1; i <= decision_diagram->topLevel(); i++) {
        for (auto& it : table[i]) {
            it.all[it.get_job()->job] = 1;
            it.all |= table.node(it.branch[1]).all;
            it.all &= table.node(it.branch[0]).all;
        }
    }

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (table.node(it.branch[1]).all[it.get_job()->job]) {
                removed_edges = true;
                it.calc_yes = false;
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        std::cout << "removing edges based on bottum-up iteration\n";
        std::cout << "Number edges removed bottum-up iteration = "
                  << nb_edges_removed_tmp << "\n";
        std::cout << "Number edges removed total = " << nb_removed_edges
                  << "\n";
        remove_layers();
        remove_edges();
        init_table();
        cleanup_arcs();
        // continue;
    }
}

void PricerSolverBdd::check_infeasible_arcs() {
    /** init table */
    bool               removed_edges = false;
    int                nb_edges_removed_tmp = 0;
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    for (int i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{nb_jobs, 0};
            it.calc_yes = true;
        }
    }

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto& n0 = table.node(it.branch[0]);
            n0.all &= it.all;
            auto& n1 = table.node(it.branch[1]);
            n1.all[it.get_job()->job] = 1;
        }
    }

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (!it.all.empty() &&
                it.all.find_first() != boost::dynamic_bitset<>::npos) {
                auto index = it.all.find_first();

                int  max = value_diff_Fij(it.get_weight(), it.get_job(),
                                         (Job*)g_ptr_array_index(jobs, index));
                bool index_bool = (index < (size_t)it.get_job()->job);
                while (index != boost::dynamic_bitset<>::npos && max <= 0) {
                    index = it.all.find_next(index);
                    index_bool = (index < (size_t)it.get_job()->job);
                    if (index != boost::dynamic_bitset<>::npos) {
                        int a = value_diff_Fij(
                            it.get_weight(), it.get_job(),
                            (Job*)g_ptr_array_index(jobs, index));
                        if (a > max) {
                            max = a;
                        }
                    }
                }

                if (max < 0 || (max == 0 && index_bool)) {
                    removed_edges = true;
                    it.calc_yes = false;
                    nb_removed_edges++;
                    nb_edges_removed_tmp++;
                }
            }
        }
    }

    if (removed_edges) {
        std::cout << "removing edges based on order\n";
        std::cout << "Number edges removed order = " << nb_edges_removed_tmp
                  << "\n";
        std::cout << "Number edges removed total = " << nb_removed_edges
                  << "\n";
        remove_layers();
        remove_edges();
        init_table();
        cleanup_arcs();
    }
}

void PricerSolverBdd::equivalent_paths_filtering() {
    /** init table */
    bool               removed_edges = false;
    int                nb_edges_removed_tmp = 0;
    EdgeTypeAccessor   edge_type_list(get(boost::edge_weight_t(), mip_graph));
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{nb_jobs, 0};
            it.calc_yes = true;
            it.calc_no = true;
            auto& n0 = table.node(it.branch[0]);
            n0.in_degree_0++;
            auto& n1 = table.node(it.branch[1]);
            n1.in_degree_1++;
        }
    }

    std::vector<int> vertices;

    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (it.in_degree_1 + it.in_degree_0 >= 2) {
                vertices.push_back(it.key);
            }
        }
    }

    for (auto& it : vertices) {
        NodeId       start_v = get(boost::vertex_name_t(), mip_graph, it);
        list<NodeId> queue;
        size_t       num_vertices = boost::num_vertices(mip_graph);

        queue.push_back(start_v);
        std::unique_ptr<bool[]> visited(new bool[num_vertices]());
        std::unique_ptr<bool[]> edge_visited(new bool[num_vertices]());
        std::unique_ptr<boost::dynamic_bitset<>[]> all(
            new dynamic_bitset<>[num_vertices]);
        std::unique_ptr<int[]> C(new int[num_vertices]);
        for (size_t i = 0; i < num_vertices; i++) {
            all[i] = dynamic_bitset<>(nb_jobs, 0);
            C[i] = 0;
        }

        auto& tmp_n = table.node(start_v);
        visited[tmp_n.key];
        auto stop = false;

        while (!queue.empty()) {
            auto currVertex = queue.front();
            queue.pop_front();
            auto it = boost::in_edges(table.node(currVertex).key, mip_graph);

            for (; it.first != it.second; it.first++) {
                NodeId adjVertex = get(boost::vertex_name_t(), mip_graph,
                                       source(*it.first, mip_graph));
                auto&  n = table.node(adjVertex);
                auto   high = edge_type_list[*it.first];

                if (!visited[n.key]) {
                    visited[n.key] = true;
                    queue.push_back(adjVertex);
                    if (high) {
                        auto& tmp_node = table.node(n.branch[1]);
                        all[n.key] |= all[tmp_node.key];
                        all[n.key][n.get_job()->job] = 1;
                        C[n.key] = C[tmp_node.key] +
                                   value_Fj(tmp_node.get_weight(), n.get_job());
                        edge_visited[n.key] = true;
                    } else {
                        auto& tmp_node = table.node(n.branch[0]);
                        all[n.key] |= all[tmp_node.key];
                        C[n.key] = C[tmp_node.key];
                    }
                } else {
                    dynamic_bitset<> tmp;
                    int              tmp_C;
                    if (high) {
                        auto& tmp_node = table.node(n.branch[1]);
                        tmp = all[tmp_node.key];
                        tmp[n.get_job()->job] = 1;
                        tmp_C = C[tmp_node.key] +
                                value_Fj(tmp_node.get_weight(), n.get_job());
                    } else {
                        auto& tmp_node = table.node(n.branch[0]);
                        tmp = all[tmp_node.key];
                        tmp_C = C[tmp_node.key];
                    }

                    if (all[n.key] == tmp) {
                        NodeId cur;
                        NodeId prev = adjVertex;
                        if (high) {
                            if (tmp_C > C[n.key]) {
                                cur = n.branch[1];
                            } else {
                                cur = n.branch[0];
                            }
                        } else {
                            if (tmp_C > C[n.key]) {
                                cur = n.branch[0];
                            } else {
                                cur = n.branch[1];
                            }
                        }

                        while (cur != start_v) {
                            auto& node = table.node(cur);
                            if (node.in_degree_1 + node.in_degree_0 > 1) {
                                break;
                            }
                            prev = cur;
                            assert(cur != NodeId(0, 0));
                            assert(cur != NodeId(0, 1));
                            if (edge_visited[node.key]) {
                                cur = node.branch[1];
                            } else {
                                cur = node.branch[0];
                            }
                        }

                        auto& node_delete = table.node(prev);
                        if (edge_visited[node_delete.key] && cur == start_v) {
                            node_delete.calc_yes = false;
                            removed_edges = true;
                            nb_edges_removed_tmp++;
                        } else if (cur == start_v) {
                            node_delete.calc_no = false;
                            removed_edges = true;
                            nb_edges_removed_tmp++;
                        }
                    }
                    stop = true;
                    break;
                }
            }

            if (stop) {
                break;
            }
        }
    }

    if (removed_edges) {
        std::cout << "Number of edges removed by equivalent_path_filtering = "
                  << nb_edges_removed_tmp << "\n"
                  << "Number of edges removed in total = "
                  << "\n";

        remove_layers();
        remove_edges();
        cleanup_arcs();
        init_table();
        construct_mipgraph();
    }
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
                                                int              num_columns) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    std::fill(lp_x.get(), lp_x.get() + get_nb_edges(), 0);
    for (int i = 0; i < num_columns; ++i) {
        if (columns[i] > 0.00001) {
            size_t       counter = 0;
            ScheduleSet* tmp =
                (ScheduleSet*)g_ptr_array_index(schedule_sets, i);
            NodeId tmp_nodeid(decision_diagram->root());

            while (tmp_nodeid > 1) {
                Job* tmp_j = nullptr;

                if (counter < tmp->job_list->len) {
                    tmp_j = (Job*)g_ptr_array_index(tmp->job_list, counter);
                }

                NodeBdd<>& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    lp_x[tmp_node.high_edge_key] += columns[i];
                    tmp_nodeid = tmp_node.branch[1];
                    counter++;
                } else {
                    lp_x[tmp_node.low_edge_key] += columns[i];
                    tmp_nodeid = tmp_node.branch[0];
                }
            }
        }
    }

    ColorWriterEdgeX  edge_writer(mip_graph, lp_x.get());
    ColorWriterVertex vertex_writer(mip_graph, table);
    string            file_name = "lp_solution_" + problem_name + "_" +
                       std::to_string(num_machines) + ".gv";
    std::ofstream outf(file_name);
    boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    outf.close();
}

void PricerSolverBdd::project_solution(Solution* sol) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    // double*            x = new double[num_edges(mip_graph)]{};
    std::fill(solution_x.get(), solution_x.get() + get_nb_edges(), 0.0);

    for (int i = 0; i < sol->nb_machines; ++i) {
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
                solution_x[tmp_node.high_edge_key] += 1.0;
                tmp_nodeid = tmp_node.branch[1];
                counter++;
            } else {
                solution_x[tmp_node.low_edge_key] += 1.0;
                tmp_nodeid = tmp_node.branch[0];
            }
        }
    }
}

void PricerSolverBdd::represent_solution(Solution* sol) {
    project_solution(sol);
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    ColorWriterEdgeX   edge_writer(mip_graph, solution_x.get());
    ColorWriterVertex  vertex_writer(mip_graph, table);
    string             file_name =
        "solution_" + problem_name + "_" + std::to_string(num_machines) + ".gv";
    std::ofstream outf(file_name);
    boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    outf.close();
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

void PricerSolverBdd::make_schedule_set_feasible(GPtrArray* set) {}

void PricerSolverBdd::disjunctive_inequality(double* x, Solution* sol) {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();
    int                branch_key = -1;
    std::unique_ptr<GRBModel> model_inequality(new GRBModel(*env));
    EdgeTypeAccessor  edge_type_list(get(boost::edge_weight_t(), mip_graph));
    EdgeVarAccessor   edge_var_list(get(boost::edge_weight2_t(), mip_graph));
    EdgeIndexAccessor edge_index_list(get(boost::edge_index_t(), mip_graph));
    VarsNodeAccessor  node_var_list(get(boost::vertex_distance_t(), mip_graph));
    NodeIdAccessor    node_id_list(get(boost::vertex_name_t(), mip_graph));
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
        GRBVar s =
            model_inequality->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        GRBVar t =
            model_inequality->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            edge_var_list[*it.first].alpha = model_inequality->addVar(
                -GRB_INFINITY, GRB_INFINITY,
                solution_x[edge_index_list[*it.first]], GRB_CONTINUOUS);
        }

        for (auto it = vertices(mip_graph); it.first != it.second; it.first++) {
            node_var_list[*it.first].omega[0] = model_inequality->addVar(
                -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            node_var_list[*it.first].omega[1] = model_inequality->addVar(
                -GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        std::unique_ptr<GRBVar[]> pi_0(new GRBVar[nb_jobs]);
        std::unique_ptr<GRBVar[]> pi_1(new GRBVar[nb_jobs]);

        for (int j = 0; j < nb_jobs; j++) {
            pi_0[j] = model_inequality->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0,
                                               GRB_CONTINUOUS);
            pi_1[j] = model_inequality->addVar(-GRB_INFINITY, GRB_INFINITY, 0.0,
                                               GRB_CONTINUOUS);
        }

        GRBVar alpha = model_inequality->addVar(-GRB_INFINITY, GRB_INFINITY,
                                                -1.0, GRB_CONTINUOUS);
        model_inequality->update();
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

        for (int j = 0; j < nb_jobs; j++) {
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
            model_inequality->addConstrs(constraints_0.get(), sense.get(),
                                         rhs.get(), nullptr, num_edges + 1));
        std::unique_ptr<GRBConstr[]> constrs1(
            model_inequality->addConstrs(constraints_1.get(), sense.get(),
                                         rhs.get(), nullptr, num_edges + 1));
        model_inequality->addConstr(normalization, GRB_EQUAL, -1.0);
        model_inequality->update();
        model_inequality->optimize();
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
            std::cout << nb_layers - *i << " ";
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

size_t PricerSolverBdd::get_nb_edges() {
    return num_edges(mip_graph);
}

size_t PricerSolverBdd::get_nb_vertices() {
    return num_vertices(mip_graph);
}

int PricerSolverBdd::get_num_layers() {
    return decision_diagram->topLevel();
}

void PricerSolverBdd::print_num_paths() {}
