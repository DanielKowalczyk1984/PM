#include "PricerSolverZdd.hpp"
#include <fmt/core.h>
#include <NodeBddStructure.hpp>
#include <cstddef>
#include "OptimalSolution.hpp"
#include "PricerConstruct.hpp"
#include "boost/graph/graphviz.hpp"
#include "gurobi_c.h"

PricerSolverZdd::PricerSolverZdd(GPtrArray*  _jobs,
                                 int         _num_machines,
                                 GPtrArray*  _ordered_jobs,
                                 const char* _p_name,
                                 double      _ub)
    : PricerSolverBase(_jobs, _num_machines, _p_name, _ub),
      ordered_jobs(_ordered_jobs)

{
    /**
     * Construction of decision diagram
     */
    PricerConstruct ps(ordered_jobs);
    decision_diagram = std::make_unique<DdStructure<NodeZdd<>>>(ps);
    remove_layers_init();
    // decision_diagram->compressBdd();
    // decision_diagram->reduceZdd();
    size_graph = decision_diagram->size();
    init_table();
    construct_mipgraph();
    lp_x = std::vector<double>(num_edges(mip_graph), 0.0);
    solution_x = std::vector<double>(num_edges(mip_graph), 0.0);
}

void PricerSolverZdd::construct_mipgraph() {
    mip_graph.clear();
    auto&             table = *(decision_diagram->getDiagram());
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
        auto job = (static_cast<job_interval_pair*>(
                        g_ptr_array_index(ordered_jobs, ordered_jobs->len - i)))
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
    auto& table = *(decision_diagram->getDiagram());
    /** init table */
    auto&     n = table.node(decision_diagram->root());
    std::span span_ordered_job{ordered_jobs->pdata, ordered_jobs->len};
    n.add_sub_node(0, decision_diagram->root(), true, false);
    n.set_node_id(decision_diagram->root());
    if (table.node(decision_diagram->root()).list.size() > 1) {
        table.node(decision_diagram->root()).list.pop_back();
    }

    for (size_t i = decision_diagram->topLevel(); i >= 0; i--) {
        size_t layer = ordered_jobs->len - i;
        auto*  tmp_pair =
            static_cast<job_interval_pair*>(span_ordered_job[layer]);

        for (auto& it : table[i]) {
            if (i != 0) {
                it.set_job(tmp_pair->j);
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
                it.set_job(nullptr);
            }
        }
    }
}

OptimalSolution<double> PricerSolverZdd::farkas_pricing(
    [[maybe_unused]] double* pi) {
    OptimalSolution<double> sol;

    return sol;
}

void PricerSolverZdd::remove_layers_init() {
    int   first_del = -1;
    int   last_del = -1;
    int   it = 0;
    auto& table = *(decision_diagram->getDiagram());

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

    fmt::print("The new number of layers = {}\n", ordered_jobs->len);
}

void PricerSolverZdd::remove_layers() {
    int   first_del = -1;
    int   last_del = -1;
    int   it = 0;
    auto& table = *(decision_diagram->getDiagram());

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        bool remove_layer = true;

        for (auto& iter : table[i]) {
            auto end =
                std::remove_if(iter.list.begin(), iter.list.end(),
                               [](std::shared_ptr<SubNodeZdd<>> const& n) {
                                   return !(n->calc_yes);
                               });
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

    fmt::print("The new number of layers = {}\n", ordered_jobs->len);
}

void PricerSolverZdd::remove_edges() {
    // decision_diagram->compressBdd();
    // decision_diagram->reduceZdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram->size();
    fmt::print("The new size of ZDD = {}\n", size_graph);
}

void PricerSolverZdd::build_mip() {
    try {
        fmt::print("Building Mip model for the extended formulation:\n");
        NodeIdAccessor vertex_nodeid_list(
            get(boost::vertex_name_t(), mip_graph));
        NodeMipIdAccessor vertex_mip_id_list(
            get(boost::vertex_degree_t(), mip_graph));
        EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(), mip_graph));
        EdgeVarAccessor  edge_var_list(get(boost::edge_weight2_t(), mip_graph));

        /** Constructing variables */
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            if (edge_type_list[*it.first]) {
                auto& n = get(boost::vertex_color_t(), mip_graph,
                              source(*it.first, mip_graph));
                Job*  job = n->get_job();

                double cost = value_Fj(n->weight + job->processing_time, job);
                edge_var_list[*it.first].x =
                    model.addVar(0.0, 1.0, cost, GRB_CONTINUOUS);
            } else {
                edge_var_list[*it.first].x = model.addVar(
                    0.0, static_cast<double>(convex_rhs), 0.0, GRB_CONTINUOUS);
            }
        }

        model.update();
        /** Assignment constraints */
        std::vector<GRBLinExpr> assignment(convex_constr_id, GRBLinExpr());
        std::vector<char>       sense(convex_constr_id, GRB_GREATER_EQUAL);
        std::vector<double>     rhs(convex_constr_id, 1.0);

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto high = edge_type_list[*it.first];

            if (high) {
                auto& n = get(boost::vertex_color_t(), mip_graph,
                              source(*it.first, mip_graph));
                assignment[n->get_job()->job] += edge_var_list[*it.first].x;
            }
        }

        std::unique_ptr<GRBConstr> assignment_constrs(
            model.addConstrs(assignment.data(), sense.data(), rhs.data(),
                             nullptr, convex_constr_id));
        model.update();
        /** Flow constraints */
        size_t num_vertices = boost::num_vertices(mip_graph);
        // boost::num_vertices(mip_graph) - table[0][1].list.size() + 1;
        std::vector<GRBLinExpr> flow_conservation_constr(num_vertices,
                                                         GRBLinExpr());
        std::vector<char>       sense_flow(num_vertices, GRB_EQUAL);
        std::vector<double>     rhs_flow(num_vertices, 0.0);

        for (auto it = vertices(mip_graph); it.first != it.second; ++it.first) {
            const auto node_id = vertex_nodeid_list[*it.first];
            const auto vertex_key = vertex_mip_id_list[*it.first];
            auto       out_edges_it = boost::out_edges(*it.first, mip_graph);

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
                rhs_flow[vertex_key] = static_cast<double>(-convex_rhs);
            } else if (node_id.row() == 0) {
                rhs_flow[vertex_key] = static_cast<double>(convex_rhs);
            }
        }

        std::unique_ptr<GRBConstr> flow_constrs(
            model.addConstrs(flow_conservation_constr.data(), sense_flow.data(),
                             rhs_flow.data(), nullptr, num_vertices));
        model.update();
        model.write("zdd_" + problem_name + "_" + std::to_string(convex_rhs) +
                    ".lp");
        model.optimize();
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
    // decision_diagram->compressBdd();
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
    auto&     table = *(decision_diagram->getDiagram());
    std::span aux_cols{columns, schedule_sets->len};
    std::span aux_sets{schedule_sets->pdata, schedule_sets->len};
    std::fill(lp_x.begin(), lp_x.end(), 0.0);
    for (int i = 0; i < num_columns; ++i) {
        if (aux_cols[i] > EPS_SOLVER) {
            size_t    counter = 0;
            auto*     tmp = static_cast<ScheduleSet*>(aux_sets[i]);
            std::span aux_jobs{tmp->job_list->pdata, tmp->job_list->pdata};
            NodeId    tmp_nodeid(decision_diagram->root());
            std::shared_ptr<SubNodeZdd<>> tmp_sub_node =
                table.node(tmp_nodeid).list[0];

            while (tmp_nodeid > 1) {
                Job* tmp_j{};

                if (counter < tmp->job_list->len) {
                    tmp_j = static_cast<Job*>(aux_jobs[counter]);
                } else {
                    tmp_j = nullptr;
                }

                NodeZdd<>& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    lp_x[tmp_sub_node->high_edge_key] += aux_cols[i];
                    tmp_nodeid = tmp_node.branch[1];
                    tmp_sub_node = tmp_sub_node->y;
                    counter++;
                } else {
                    lp_x[tmp_sub_node->low_edge_key] += aux_cols[i];
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

bool PricerSolverZdd::check_schedule_set(GPtrArray* set) {
    guint     weight = 0;
    std::span aux_jobs{set->pdata, set->len};
    auto&     table = *(decision_diagram->getDiagram());
    NodeId    tmp_nodeid(decision_diagram->root());

    for (unsigned j = 0; j < set->len; ++j) {
        Job* tmp_j = static_cast<Job*>(aux_jobs[j]);
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

void PricerSolverZdd::make_schedule_set_feasible(
    [[maybe_unused]] GPtrArray* set) {}

void PricerSolverZdd::iterate_zdd() {
    DdStructure<NodeZdd<double>>::const_iterator it = decision_diagram->begin();

    for (; it != decision_diagram->end(); ++it) {
        auto i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            std::cout << ordered_jobs->len - *i << " ";
        }

        std::cout << '\n';
    }
}

void PricerSolverZdd::create_dot_zdd(const char* name) {
    std::ofstream file;
    file.open(name);
    // decision_diagram->dumpDot(file);
    file.close();
}

void PricerSolverZdd::print_number_nodes_edges() {
    fmt::print("removed edges = %d, removed nodes = {}\n", nb_removed_edges,
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
