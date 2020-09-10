#include "PricerSolverBdd.hpp"
#include <fmt/core.h>
#include <algorithm>
#include <boost/concept_archetype.hpp>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iostream>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>
#include "MipGraph.hpp"
#include "ModelInterface.hpp"
#include "NodeBdd.hpp"
#include "NodeBddStructure.hpp"
#include "NodeId.hpp"
#include "OptimalSolution.hpp"
#include "PricerConstruct.hpp"
#include "ZeroHalfCuts.hpp"
#include "boost/graph/graphviz.hpp"
#include "gurobi_c++.h"
#include "gurobi_c.h"
#include "interval.h"
#include "job.h"
#include "lp.h"
#include "scheduleset.h"
#include "util.h"
#include "wctprivate.h"

using namespace std;

PricerSolverBdd::PricerSolverBdd(GPtrArray*  _jobs,
                                 int         _num_machines,
                                 GPtrArray*  _ordered_jobs,
                                 const char* p_name,
                                 int         _Hmax,
                                 int*        _take_jobs,
                                 double      _UB)
    : PricerSolverBase(_jobs, _num_machines, p_name, _UB),
      size_graph(0),
      nb_removed_edges(0),
      nb_removed_nodes(0),
      ordered_jobs(_ordered_jobs),
      nb_layers(ordered_jobs->len),
      node_ids(convex_constr_id, vector<std::weak_ptr<NodeId>>(_Hmax + 1)),
      original_model(reformulation_model),
      H_min(0),
      H_max(_Hmax)

{
    /**
     * Construction of decision diagram
     */
    if (_take_jobs) {
        PricerConstructTI ps(ordered_jobs, _take_jobs, _Hmax);
        decision_diagram = std::make_unique<DdStructure<>>(ps);
    } else {
        PricerConstruct ps(ordered_jobs);
        decision_diagram = std::make_unique<DdStructure<>>(ps);
    }
    remove_layers_init();
    decision_diagram->compressBdd();
    size_graph = decision_diagram->size();
    init_table();
    calculate_H_min();
    cleanup_arcs();
    // check_infeasible_arcs();
    bottum_up_filtering();
    topdown_filtering();
    construct_mipgraph();
    init_coeff_constraints();
    std::cout << "Ending construction\n";
    solution_x = std::unique_ptr<double[]>(new double[get_nb_edges()]);
}

void PricerSolverBdd::calculate_H_min() {
    auto p_sum = 0.0;
    auto duration = g_ptr_array_new();
    for (auto j = 0; j < convex_constr_id; j++) {
        auto job = (Job*)g_ptr_array_index(jobs, j);
        g_ptr_array_add(duration, job);
        p_sum += job->processing_time;
    }
    g_ptr_array_sort(duration, g_compare_duration);

    auto m = 0;
    auto tmp = p_sum;
    auto i = convex_constr_id;
    do {
        auto job = (Job*)g_ptr_array_index(duration, i - 1);
        tmp -= job->processing_time;
        m++;
        i--;
    } while (m < convex_rhs - 1);

    H_min = static_cast<int>(floor(tmp / convex_rhs));

    g_ptr_array_free(duration, TRUE);
}

void PricerSolverBdd::construct_mipgraph() {
    mip_graph.clear();
    auto& table = decision_diagram->getDiagram().privateEntity();
    auto  vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    auto  edge_type_list(get(boost::edge_weight_t(), mip_graph));

    for (auto i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto j = 0u; j < table[i].size(); j++) {
            if (NodeId(i, j) != 0
                // && (table[i][j].calc_yes || table[i][j].calc_no)
            ) {
                table[i][j].key = add_vertex(mip_graph);
                vertex_nodeid_list[table[i][j].key] = NodeId(i, j);
            }
        }
    }

    auto count = 0;

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
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

    // std::cout << "Number of vertices = " << num_vertices(mip_graph) << '\n';
    // std::cout << "Number of edges = " << num_edges(mip_graph) << '\n';
}

void PricerSolverBdd::init_coeff_constraints() {
    auto& table = decision_diagram->getDiagram().privateEntity();

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (auto c = 0; c < reformulation_model.get_nb_constraints();
                 c++) {
                if (c == convex_constr_id) {
                    continue;
                }

                auto            constr = reformulation_model.get_constraint(c);
                VariableKeyBase key_aux(it.get_nb_job(), it.get_weight());
                auto            coeff = constr->get_var_coeff(&key_aux);
                if (fabs(coeff) > 1e-10) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff, 0.0, c)};
                    original_model.add_coeff_list(c, ptr_coeff);
                    it.add_coeff_list(ptr_coeff, 1);
                }
            }
        }
    }

    auto&           root = decision_diagram->root();
    auto&           root_node = table.node(root);
    VariableKeyBase key_aux(root_node.get_nb_job(), root_node.get_weight(),
                            true);
    auto            constr = original_model.get_constraint(convex_constr_id);
    auto            coeff = constr->get_var_coeff(&key_aux);
    if (fabs(coeff) > 1e-10) {
        auto ptr_coeff_high{std::make_shared<BddCoeff>(
            root_node.get_nb_job(), root_node.get_weight(), coeff, true, true)};
        original_model.add_coeff_list(convex_constr_id, ptr_coeff_high);
        auto ptr_coeff_low{std::make_shared<BddCoeff>(root_node.get_nb_job(),
                                                      root_node.get_weight(),
                                                      coeff, false, true)};
        original_model.add_coeff_list(convex_constr_id, ptr_coeff_low);
    }
}

void PricerSolverBdd::update_coeff_constraints() {
    int  nb_constr = original_model.get_nb_constraints();
    auto nb_new_constr = reformulation_model.get_nb_constraints() - nb_constr;

    for (auto j = nb_constr; j < reformulation_model.get_nb_constraints();
         j++) {
        original_model.add_constraint(reformulation_model.get_constraint(j));
    }

    auto& table = decision_diagram->getDiagram().privateEntity();
    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (int c = 0; c < nb_new_constr; c++) {
                auto     constr = original_model.get_constraint(nb_constr + c);
                BddCoeff key_high{it.get_nb_job(), it.get_weight(), 0.0, 0.0,
                                  nb_constr + c};
                auto     coeff_high = constr->get_var_coeff(&key_high);
                if (fabs(coeff_high) > 1e-10) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff_high, 0.0,
                        nb_constr + c)};
                    original_model.add_coeff_list(c, ptr_coeff);
                    it.add_coeff_list(ptr_coeff, 1);
                }

                BddCoeff key_low{it.get_nb_job(),
                                 it.get_weight(),
                                 0.0,
                                 0.0,
                                 nb_constr + c,
                                 false};
                auto     coeff_low = constr->get_var_coeff(&key_high);
                if (fabs(coeff_low) > 1e-10) {
                    std::shared_ptr<BddCoeff> ptr_coeff{
                        std::make_shared<BddCoeff>(it.get_nb_job(),
                                                   it.get_weight(), coeff_low,
                                                   0.0, nb_constr + c, false)};
                    original_model.add_coeff_list(c, ptr_coeff);
                    it.add_coeff_list(ptr_coeff, 0);
                }
            }
        }
    }
}

void PricerSolverBdd::init_table() {
    auto& table = decision_diagram->getDiagram().privateEntity();
    /** init table */
    auto& root = table.node(decision_diagram->root());
    root.init_node(0, true);
    root.all = boost::dynamic_bitset<>{convex_constr_id, 0};

    for (auto i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto it = 0u; it < table[i].size(); it++) {
            if (i != 0) {
                auto layer = nb_layers - i;
                auto tmp_pair = reinterpret_cast<job_interval_pair*>(
                    g_ptr_array_index(ordered_jobs, layer));
                auto& node = table.node(NodeId(i, it));
                auto  aux_job = tmp_pair->j;
                auto  w = node.get_weight();
                auto  p = aux_job->processing_time;

                auto& n0 = table.node(node.branch[0]);
                auto& n1 = table.node(node.branch[1]);

                node.ptr_node_id = std::make_shared<NodeId>(i, it);
                node.set_job(aux_job);
                node.child[0] = n0.init_node(w);
                node.child[1] = n1.init_node(w + p);
                node.cost[0] = 0.0;
                node.cost[1] = value_Fj(w + p, aux_job);

                node_ids[aux_job->job][w] = node.ptr_node_id;

                n0.in_degree_0++;
                n0.in_edges[0].push_back(node.ptr_node_id);
                n1.in_degree_1++;
                n1.in_edges[1].push_back(node.ptr_node_id);

                auto iter = t_out.find(w);
                if (iter == t_out.end()) {
                    t_out[w] =
                        std::vector<std::weak_ptr<NodeId>>{node.ptr_node_id};
                } else {
                    iter->second.push_back(node.ptr_node_id);
                }

                iter = t_in.find(w + p);
                if (iter == t_in.end()) {
                    t_in[w + p] =
                        std::vector<std::weak_ptr<NodeId>>{node.ptr_node_id};
                } else {
                    iter->second.push_back(node.ptr_node_id);
                }

            } else {
                auto& node = table.node(NodeId(i, it));
                node.set_job(nullptr);
            }
        }
    }
}

void PricerSolverBdd::insert_constraints_lp(NodeData* pd) {
    wctlp_get_nb_rows(pd->RMP, &(pd->nb_rows));
    int nb_new_constraints =
        reformulation_model.get_nb_constraints() - pd->nb_rows;

    assert(nb_new_constraints <=
           (pd->id_pseudo_schedules - pd->id_next_var_cuts));
    std::vector<int>    starts(nb_new_constraints + 1);
    std::vector<char>   sense(nb_new_constraints);
    std::vector<double> rhs(nb_new_constraints);
    std::vector<int>    column_ind;
    std::vector<double> coeff;

    int pos = 0;
    for (int c = 0; c < nb_new_constraints; c++) {
        auto constr = reformulation_model.get_constraint(pd->nb_rows + c);

        sense[c] = constr->get_sense();
        starts[c] = pos;
        rhs[c] = constr->get_rhs();

        if (rhs[c] != 0.0) {
            pos++;
            column_ind.push_back(pd->id_next_var_cuts++);
            if (sense[c] == '>') {
                coeff.push_back(1.0);
            } else {
                coeff.push_back(-1.0);
            }
        }

        for (auto i = 0u; i < pd->localColPool->len; i++) {
            auto aux_schedule_set =
                (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
            auto  jobs_list = aux_schedule_set->job_list;
            auto& table = decision_diagram->getDiagram().privateEntity();
            auto  tmp_nodeid(decision_diagram->root());

            auto coeff_val = 0.0;
            auto counter = 0u;
            while (tmp_nodeid > 1) {
                auto& tmp_node = table.node(tmp_nodeid);
                Job*  tmp_j = nullptr;

                if (counter < jobs_list->len) {
                    tmp_j = (Job*)g_ptr_array_index(jobs_list, counter);
                }

                VariableKeyBase key(tmp_node.get_nb_job(),
                                    tmp_node.get_weight(),
                                    tmp_j == tmp_node.get_job());
                if (key.get_high()) {
                    coeff_val += constr->get_var_coeff(&key);
                    tmp_nodeid = tmp_node.branch[1];
                    counter++;
                } else {
                    coeff_val += constr->get_var_coeff(&key);
                    tmp_nodeid = tmp_node.branch[0];
                }
            }

            assert(tmp_nodeid == 1);

            if (fabs(coeff_val) > 1e-6) {
                column_ind.push_back(pd->id_pseudo_schedules + i);
                coeff.push_back(coeff_val);
                pos++;
            }
        }
    }

    starts[nb_new_constraints] = pos;

    wctlp_addrows(pd->RMP, nb_new_constraints, coeff.size(), starts.data(),
                  column_ind.data(), coeff.data(), sense.data(), rhs.data(),
                  nullptr);
    wctlp_get_nb_rows(pd->RMP, &(pd->nb_rows));

    vector<double> new_values(nb_new_constraints, 0.0);
    vector<int>    new_values_int(nb_new_constraints, 0);
    g_array_append_vals(pd->pi, new_values.data(), new_values.size());
    g_array_append_vals(pd->rhs, new_values.data(), new_values.size());
    wctlp_get_rhs(pd->RMP, &g_array_index(pd->rhs, double, 0));
    g_array_append_vals(pd->lhs_coeff, new_values.data(), new_values.size());
    g_array_append_vals(pd->id_row, new_values_int.data(),
                        new_values_int.size());
    g_array_append_vals(pd->coeff_row, new_values.data(), new_values.size());
}

double PricerSolverBdd::compute_reduced_cost(const OptimalSolution<>& sol,
                                             double*                  pi,
                                             double*                  lhs) {
    double result = sol.cost;
    auto&  table = *decision_diagram->getDiagram();
    auto   tmp_nodeid(decision_diagram->root());
    auto   counter = 0u;

    std::fill(lhs, lhs + reformulation_model.get_nb_constraints(), 0.0);
    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        Job*  tmp_j = nullptr;

        if (counter < sol.jobs->len) {
            tmp_j = (Job*)g_ptr_array_index(sol.jobs, counter);
        }

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node.branch[1];
            counter++;
            auto dual = pi[key.get_j()];
            auto constr = reformulation_model.get_constraint(key.get_j());
            auto coeff = constr->get_var_coeff(&key);

            if (fabs(coeff) > 1e-10) {
                result -= coeff * dual;
                lhs[key.get_j()] += coeff;
            }
        } else {
            tmp_nodeid = tmp_node.branch[0];
        }

        for (int c = convex_constr_id + 1;
             c < reformulation_model.get_nb_constraints(); c++) {
            if (c == convex_constr_id) {
                continue;
            }
            auto dual = pi[c];
            auto constr = reformulation_model.get_constraint(c);
            auto coeff = constr->get_var_coeff(&key);

            if (fabs(coeff) > 1e-10) {
                result -= coeff * dual;
                lhs[c] += coeff;
            }
        }
    }

    auto dual = pi[convex_constr_id];
    auto constr = reformulation_model.get_constraint(convex_constr_id);
    VariableKeyBase k(0, 0, true);
    auto            coeff = constr->get_var_coeff(&k);
    result -= coeff * dual;
    lhs[convex_constr_id] += coeff;

    return result;
}

double PricerSolverBdd::compute_subgradient(const OptimalSolution<>& sol,
                                            double* sub_gradient) {
    double result = sol.cost;
    auto&  table = *decision_diagram->getDiagram();
    auto   tmp_nodeid(decision_diagram->root());
    auto   counter = 0u;
    auto   nb_constraints = reformulation_model.get_nb_constraints();
    auto   convex_rhs =
        -reformulation_model.get_constraint(convex_constr_id)->get_rhs();

    for (size_t i = 0; i < nb_constraints; i++) {
        auto constr = reformulation_model.get_constraint(i);
        sub_gradient[i] = constr->get_rhs();
    }

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        Job*  tmp_j = nullptr;

        if (counter < sol.jobs->len) {
            tmp_j = (Job*)g_ptr_array_index(sol.jobs, counter);
        }

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node.branch[1];
            counter++;
            auto constr = reformulation_model.get_constraint(key.get_j());
            auto coeff = constr->get_var_coeff(&key);

            if (fabs(coeff) > 1e-10) {
                sub_gradient[key.get_j()] -= coeff * convex_rhs;
            }
        } else {
            tmp_nodeid = tmp_node.branch[0];
        }

        for (int c = convex_constr_id + 1;
             c < reformulation_model.get_nb_constraints(); c++) {
            // auto dual = pi[c];
            auto constr = reformulation_model.get_constraint(c);
            auto coeff = constr->get_var_coeff(&key);

            if (fabs(coeff) > 1e-10) {
                sub_gradient[c] -= coeff * convex_rhs;
            }
        }
    }

    sub_gradient[convex_constr_id] += convex_rhs;
    assert(sub_gradient[convex_constr_id] == 0.0);
    // sub_gradient[nb_jobs] = 0.0;
    assert(tmp_nodeid == 1);

    return result;
}

double PricerSolverBdd::compute_lagrange(const OptimalSolution<>& sol,
                                         double*                  pi) {
    double result = sol.cost;
    auto   dual_bound = 0.0;

    auto& table = *decision_diagram->getDiagram();
    auto  tmp_nodeid(decision_diagram->root());

    auto counter = 0u;
    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        Job*  tmp_j = nullptr;

        if (counter < sol.jobs->len) {
            tmp_j = (Job*)g_ptr_array_index(sol.jobs, counter);
        }

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            auto dual = pi[key.get_j()];
            auto constr = reformulation_model.get_constraint(key.get_j());
            auto coeff = constr->get_var_coeff(&key);

            if (fabs(coeff) > 1e-10) {
                result -= coeff * dual;
            }

            counter++;
            tmp_nodeid = tmp_node.branch[1];
        } else {
            tmp_nodeid = tmp_node.branch[0];
        }

        for (int c = convex_constr_id + 1;
             c < reformulation_model.get_nb_constraints(); c++) {
            if (c == convex_constr_id) {
                continue;
            }
            auto dual = pi[c];
            auto constr = reformulation_model.get_constraint(c);
            auto coeff = constr->get_var_coeff(&key);

            if (fabs(coeff) > 1e-10) {
                result -= coeff * dual;
            }
        }
    }

    result = CC_MIN(0, result);

    for (int c = 0; c < reformulation_model.get_nb_constraints(); c++) {
        if (c == convex_constr_id) {
            continue;
        }
        auto dual = pi[c];
        auto constr = reformulation_model.get_constraint(c);
        auto rhs = constr->get_rhs();

        dual_bound += rhs * dual;
    }

    result = -reformulation_model.get_constraint(convex_constr_id)->get_rhs() *
             result;
    result = dual_bound + result;

    // assert(tmp_nodeid == 1);

    return result;
}
void PricerSolverBdd::remove_layers_init() {
    auto  first_del = -1;
    auto  last_del = -1;
    auto  it = 0;
    auto& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
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
    fmt::print("{0: <{2}}{1}\n", "The new number of layers", nb_layers, 60);
}

void PricerSolverBdd::remove_layers() {
    auto  first_del = -1;
    auto  last_del = -1;
    auto  it = 0;
    auto& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary layers of the bdd */
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        auto remove = true;

        for (auto& iter : table[i]) {
            if (iter.calc_yes) {
                remove = false;
            } else {
                auto& cur_node_1 = iter.branch[1];
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
    fmt::print("{0: <{2}}{1}\n", "The new number of layers", nb_layers, 60);
}

void PricerSolverBdd::remove_edges() {
    auto& table = decision_diagram->getDiagram().privateEntity();

    /** remove the unnecessary nodes of the bdd */
    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& iter : table[i]) {
            if (!iter.calc_yes) {
                NodeId& cur_node_1 = iter.branch[1];
                iter.ptr_node_id.reset();
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
    // printf("The new size of BDD \t\t\t\t= %lu\n", size_graph);
    // std::cout
    //     << "-------------------------------------------------------------\n";
}

void PricerSolverBdd::print_representation_file() {
    auto& table = decision_diagram->getDiagram().privateEntity();
    auto  vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    auto  edge_type_list(get(boost::edge_weight_t(), mip_graph));
    auto  edge_index_list(get(boost::edge_index_t(), mip_graph));
    auto  index_edge{std::make_unique<std::vector<int>[]>(convex_constr_id)};

    auto outfile_file_mip_str =
        problem_name + "_" + std::to_string(convex_rhs) + ".txt";
    std::ofstream out_file_mip(outfile_file_mip_str);

    out_file_mip << boost::num_vertices(mip_graph) << " "
                 << boost::num_edges(mip_graph) << " " << convex_constr_id
                 << " " << convex_rhs << "\n\n";

    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto& head = table.node(get(boost::vertex_name_t(), mip_graph,
                                    source(*it.first, mip_graph)));
        auto& n = table.node(get(boost::vertex_name_t(), mip_graph,
                                 target(*it.first, mip_graph)));
        auto  high = edge_type_list[*it.first];
        if (high) {
            double cost =
                value_Fj(head.get_weight() + head.get_job()->processing_time,
                         head.get_job());
            out_file_mip << head.key << " " << n.key << " " << cost << "\n";
            index_edge[head.get_nb_job()].push_back(edge_index_list[*it.first]);
        } else {
            out_file_mip << head.key << " " << n.key << " " << 0.0 << "\n";
        }
    }

    out_file_mip << "\n";

    for (int i = 0; i < convex_constr_id; i++) {
        out_file_mip << index_edge[i].size() << " ";
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
            out_file_mip << n.get_nb_job() << " " << n.get_weight() << "\n";
        }
    }

    out_file_mip << "99 99\n";

    out_file_mip.close();
}

void PricerSolverBdd::add_inequality(std::vector<int> v1, std::vector<int> v2) {
    GRBLinExpr        expr1;
    EdgeVarAccessor   edge_var_list(get(boost::edge_weight2_t(), mip_graph));
    EdgeIndexAccessor edge_index_list(get(boost::edge_index_t(), mip_graph));
    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        if (std::find(v1.begin(), v1.end(), edge_index_list[*it.first]) !=
            v1.end()) {
            expr1 += edge_var_list[*it.first].x;
        }

        if (std::find(v2.begin(), v2.end(), edge_index_list[*it.first]) !=
            v2.end()) {
            expr1 -= edge_var_list[*it.first].x;
        }
    }
    model->addConstr(expr1, GRB_EQUAL, 0);
}

void PricerSolverBdd::add_inequality(std::vector<int> v1) {
    GRBLinExpr expr1;
    auto       edge_var_list(get(boost::edge_weight2_t(), mip_graph));
    auto       edge_index_list(get(boost::edge_index_t(), mip_graph));
    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        if (std::find(v1.begin(), v1.end(), edge_index_list[*it.first]) !=
            v1.end()) {
            expr1 += edge_var_list[*it.first].x;
        }
    }
    model->addConstr(expr1, GRB_EQUAL, convex_rhs);
}
void PricerSolverBdd::build_mip() {
    try {
        printf("Building Mip model for the extended formulation:\n");
        auto& table = decision_diagram->getDiagram().privateEntity();
        auto  vertex_index_list(get(boost::vertex_index_t(), mip_graph));
        auto  vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
        auto  edge_type_list{get(boost::edge_weight_t(), mip_graph)};
        auto  edge_var_list{get(boost::edge_weight2_t(), mip_graph)};
        auto  edge_index_list{get(boost::edge_index_t(), mip_graph)};

        /** Constructing variables */
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            if (edge_type_list[*it.first]) {
                auto&  n = table.node(get(boost::vertex_name_t(), mip_graph,
                                         source(*it.first, mip_graph)));
                auto   C = n.get_weight() + n.get_job()->processing_time;
                double cost = value_Fj(C, n.get_job());
                edge_var_list[*it.first].x =
                    model->addVar(0.0, 1.0, cost, GRB_BINARY);
            } else {
                edge_var_list[*it.first].x = model->addVar(
                    0.0, static_cast<double>(convex_rhs), 0.0, GRB_CONTINUOUS);
            }
        }

        model->update();
        /** Assignment constraints */
        auto assignment{std::make_unique<GRBLinExpr[]>(convex_constr_id)};
        auto sense{std::make_unique<char[]>(convex_constr_id)};
        auto rhs{std::make_unique<double[]>(convex_constr_id)};

        for (unsigned i = 0; i < jobs->len; ++i) {
            sense[i] = GRB_EQUAL;
            rhs[i] = 1.0;
        }

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto high = edge_type_list[*it.first];

            if (high) {
                auto& n = table.node(get(boost::vertex_name_t(), mip_graph,
                                         source(*it.first, mip_graph)));
                assignment[n.get_nb_job()] += edge_var_list[*it.first].x;
            }
        }

        std::unique_ptr<GRBConstr[]> assignment_constrs(
            model->addConstrs(assignment.get(), sense.get(), rhs.get(), nullptr,
                              convex_constr_id));
        model->update();
        /** Flow constraints */
        auto num_vertices = boost::num_vertices(mip_graph);
        auto flow_conservation_constr{
            std::make_unique<GRBLinExpr[]>(num_vertices)};
        auto sense_flow{std::make_unique<char[]>(num_vertices)};
        auto rhs_flow(std::make_unique<double[]>(num_vertices));

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
                rhs_flow[vertex_key] = static_cast<double>(-convex_rhs);
            } else if (node_id == 1) {
                rhs_flow[vertex_key] = static_cast<double>(convex_rhs);
            } else {
                rhs_flow[vertex_key] = 0.0;
            }
        }

        std::unique_ptr<GRBConstr[]> flow_constrs(
            model->addConstrs(flow_conservation_constr.get(), sense_flow.get(),
                              rhs_flow.get(), nullptr, num_vertices));
        model->update();
        // for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        //     // edge_var_list[*it.first].x.set(GRB_DoubleAttr_PStart,
        //     // lp_x[edge_index_list[*it.first]]);
        //     edge_var_list[*it.first].x.set(
        //         GRB_DoubleAttr_Start,
        //         solution_x[edge_index_list[*it.first]]);
        // }
        model->write("original_" + problem_name + ".lp");
        auto presolve = model->presolve();
        presolve.write("presolve_" + problem_name + ".lp");
        model->optimize();

        if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            for (auto it = edges(mip_graph); it.first != it.second;
                 it.first++) {
                int index = edge_index_list[*it.first];
                solution_x[index] =
                    edge_var_list[*it.first].x.get(GRB_DoubleAttr_X);
            }

            // ColorWriterEdgeX  edge_writer(mip_graph, solution_x.get());
            // ColorWriterVertex vertex_writer(mip_graph, table);
            // string            file_name = "lp_solution_" + problem_name + "_"
            // +
            //                    std::to_string(num_machines) + ".gv";
            // std::ofstream outf(file_name);
            // boost::write_graphviz(outf, mip_graph, vertex_writer,
            // edge_writer); outf.close();
        }

        ColorWriterEdgeIndex edge_writer_index(mip_graph);
        ColorWriterVertex    vertex_writer(mip_graph, table);
        auto                 file_name =
            "index_" + problem_name + "_" + std::to_string(convex_rhs) + ".gv";
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
    bottum_up_filtering();
    topdown_filtering();
    cleanup_arcs();
    construct_mipgraph();
    // for (int c = 0; c < reformulation_model.get_nb_constraints(); c++) {
    //     if (c == nb_jobs) {
    //         continue;
    //     }

    //     auto coeff_list = original_model.get_coeff_list(c);
    //     auto iter = coeff_list->begin();

    //     while (iter != coeff_list->end()) {
    //         if (node_ids[(*iter)->get_j()][(*iter)->get_t()].lock()) {
    //             iter = coeff_list->erase(iter);
    //         } else {
    //             iter++;
    //         }
    //     }
    // }

    // auto& table = decision_diagram->getDiagram().privateEntity();
    // std::for_each(t_out.begin(), t_out.end(), [&](auto& elem) {
    //     std::cout << elem.first << ": ";
    //     for (auto& it : elem.second) {
    //         auto aux = it.lock();
    //         if (aux) {
    //             auto& node = table.node(*aux);
    //             std::cout << node.get_nb_job() << " ";
    //         }
    //     }
    //     std::cout << "\n";
    // });

    // std::cout << "----------------------------------------------------\n";
    // // NodeTableEntity<>&   table =
    // std::for_each(t_in.begin(), t_in.end(), [&](auto& elem) {
    //     std::cout << elem.first << ": ";
    //     for (auto& it : elem.second) {
    //         auto aux = it.lock();
    //         if (aux) {
    //             auto& node = table.node(*aux);
    //             std::cout << node.get_nb_job() << " ";
    //         }
    //     }
    //     std::cout << "\n";
    // });
    // decision_diagram->getDiagram().privateEntity(); ColorWriterEdgeIndex
    // edge_writer(mip_graph); ColorWriterVertex    vertex_writer(mip_graph,
    // table); string               file_name = "representation_" + problem_name
    // + "_" +
    //                    std::to_string(num_machines) + ".gv";
    // std::ofstream outf(file_name);
    // boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    // outf.close();
}

void PricerSolverBdd::cleanup_arcs() {
    NodeTableEntity<>& table = decision_diagram->getDiagram().privateEntity();

    table.node(0).backward_distance[0] = INT_MIN;
    table.node(0).backward_distance[1] = INT_MIN;
    table.node(1).backward_distance[0] = 0;
    table.node(1).backward_distance[1] = 0;
    auto removed_edges = false;
    auto nb_edges_removed_tmp = 0;

    for (auto i = 1; i <= decision_diagram->topLevel(); i++) {
        for (auto& it : table[i]) {
            it.calc_no = true;
            it.calc_yes = true;
            NodeBdd<>& cur_node_0 = table.node(it.branch[0]);
            NodeBdd<>& cur_node_1 = table.node(it.branch[1]);

            if (cur_node_0.backward_distance[0] <=
                cur_node_0.backward_distance[1]) {
                it.backward_distance[0] = cur_node_0.backward_distance[1];
            } else {
                it.backward_distance[0] = cur_node_0.backward_distance[0];
            }

            int result0 =
                cur_node_1.backward_distance[0] + it.get_job()->processing_time;
            int result1 =
                cur_node_1.backward_distance[1] + it.get_job()->processing_time;

            if (result0 <= result1) {
                it.backward_distance[1] = result1;
            } else {
                it.backward_distance[1] = result0;
            }
        }
    }
    /** remove the unnecessary nodes of the bdd */
    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& iter : table[i]) {
            if (iter.get_weight() + iter.backward_distance[0] < H_min &&
                iter.branch[0] != 0) {
                iter.calc_no = false;
                removed_edges = true;
                nb_edges_removed_tmp++;
                nb_removed_edges++;
            }

            if (iter.get_weight() + iter.backward_distance[1] < H_min) {
                iter.calc_yes = false;
                removed_edges = true;
                nb_edges_removed_tmp++;
                nb_removed_edges++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("{0: <{2}}{1}\n", "Number of edges removed by clean up",
                   nb_edges_removed_tmp, 60);
        fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                   get_nb_removed_edges(), 60);
        remove_layers();
        remove_edges();
        // init_table();
    }
}

void PricerSolverBdd::topdown_filtering() {
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = decision_diagram->getDiagram().privateEntity();
    auto& root = table.node(decision_diagram->root());
    root.init_node(0, true);
    for (auto i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
            it.calc_yes = true;
        }
    }

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto& n0 = table.node(it.branch[0]);

            if (n0.visited) {
                n0.all &= it.all;
            } else {
                n0.visited = true;
                n0.all = boost::dynamic_bitset<>{convex_constr_id, 0};
                n0.all |= it.all;
            }

            auto& n1 = table.node(it.branch[1]);

            if (n1.visited) {
                if (n1.all[it.get_nb_job()]) {
                    n1.all &= it.all;
                    n1.all[it.get_nb_job()] = 1;
                } else {
                    n1.all &= it.all;
                }
            } else {
                n1.all = boost::dynamic_bitset<>{convex_constr_id, 0};
                n1.all |= it.all;
                n1.all[it.get_nb_job()] = 1;
                n1.visited = true;
            }
        }
    }

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (it.all[it.get_nb_job()]) {
                removed_edges = true;
                it.calc_yes = false;
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        // std::cout << "Removing edges based on top-down iteration\n";
        // std::cout << "Number edges removed top-bottom \t\t= "
        //           << nb_edges_removed_tmp << "\n";
        // std::cout << "Number edges removed total \t\t\t= " <<
        // nb_removed_edges
        //           << "\n";
        remove_layers();
        remove_edges();
        cleanup_arcs();
        // init_table();
        // continue;
    }
}

void PricerSolverBdd::bottum_up_filtering() {
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = decision_diagram->getDiagram().privateEntity();
    for (int i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
            it.calc_yes = true;
        }
    }

    table.node(0).all = boost::dynamic_bitset<>{convex_constr_id, 0};
    table.node(1).all = boost::dynamic_bitset<>{convex_constr_id, 0};
    table.node(0).all.flip();

    for (auto i = 1; i <= decision_diagram->topLevel(); i++) {
        for (auto& it : table[i]) {
            it.all[it.get_nb_job()] = 1;
            it.all |= table.node(it.branch[1]).all;
            it.all &= table.node(it.branch[0]).all;
        }
    }

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (table.node(it.branch[1]).all[it.get_nb_job()]) {
                removed_edges = true;
                it.calc_yes = false;
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        // std::cout << "removing edges based on bottum-up iteration\n";
        // std::cout << "Number edges removed bottum-up iteration \t\t= "
        //           << nb_edges_removed_tmp << "\n";
        // std::cout << "Number edges removed total \t\t\t= " <<
        // nb_removed_edges
        //           << "\n";
        remove_layers();
        remove_edges();
        cleanup_arcs();
        // init_table();
        // continue;
    }
}

void PricerSolverBdd::check_infeasible_arcs() {
    /** init table */
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = decision_diagram->getDiagram().privateEntity();
    for (auto i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
            it.calc_yes = true;
            it.calc_no = true;
        }
    }

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto& n0 = table.node(it.branch[0]);
            n0.all |= it.all;
            auto& n1 = table.node(it.branch[1]);
            n1.all[it.get_nb_job()] = 1;
        }
    }

    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            if (!it.all.empty() &&
                it.all.find_first() != boost::dynamic_bitset<>::npos) {
                auto index = it.all.find_first();

                auto max = value_diff_Fij(it.get_weight(), it.get_job(),
                                          (Job*)g_ptr_array_index(jobs, index));
                // bool index_bool = (index > (size_t)it.get_nb_job());
                while (index != boost::dynamic_bitset<>::npos && max < 0) {
                    index = it.all.find_next(index);
                    if (index != boost::dynamic_bitset<>::npos) {
                        int a = value_diff_Fij(
                            it.get_weight(), it.get_job(),
                            (Job*)g_ptr_array_index(jobs, index));
                        if (a > max) {
                            max = a;
                        }
                    }
                }

                if (max < 0) {
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
        cleanup_arcs();
        // init_table();
    }
}

void PricerSolverBdd::equivalent_paths_filtering() {
    /** init table */
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto  edge_type_list(get(boost::edge_weight_t(), mip_graph));
    auto& table = decision_diagram->getDiagram().privateEntity();
    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
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
        auto         start_v = get(boost::vertex_name_t(), mip_graph, it);
        auto         num_vertices = boost::num_vertices(mip_graph);
        list<NodeId> queue;

        queue.push_back(start_v);
        std::unique_ptr<bool[]> visited(new bool[num_vertices]());
        std::unique_ptr<bool[]> edge_visited(new bool[num_vertices]());
        std::unique_ptr<boost::dynamic_bitset<>[]> all(
            new dynamic_bitset<>[num_vertices]);
        std::unique_ptr<int[]> C(new int[num_vertices]);
        for (auto i = 0u; i < num_vertices; i++) {
            all[i] = dynamic_bitset<>(convex_constr_id, 0);
            C[i] = 0;
        }

        auto& tmp_n = table.node(start_v);
        visited[tmp_n.key];
        auto stop = false;

        while (!queue.empty()) {
            auto currVertex = queue.front();
            queue.pop_front();
            auto iter = boost::in_edges(table.node(currVertex).key, mip_graph);

            for (; iter.first != iter.second; iter.first++) {
                auto  adjVertex = get(boost::vertex_name_t(), mip_graph,
                                     source(*iter.first, mip_graph));
                auto& n = table.node(adjVertex);
                auto  high = edge_type_list[*iter.first];

                if (!visited[n.key]) {
                    visited[n.key] = true;
                    queue.push_back(adjVertex);
                    if (high) {
                        auto& tmp_node = table.node(n.branch[1]);
                        all[n.key] |= all[tmp_node.key];
                        all[n.key][n.get_nb_job()] = 1;
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
                        tmp[n.get_nb_job()] = 1;
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
        // init_table();
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
    // init_table();
    cout << decision_diagram->size() << '\n';
    construct_mipgraph();
    auto&             table1 = decision_diagram->getDiagram().privateEntity();
    ColorWriterVertex vertex_writer1(mip_graph, table1);
    outf = std::ofstream("min2.gv");
    boost::write_graphviz(outf, mip_graph, vertex_writer1);
    outf.close();
}

void PricerSolverBdd::construct_lp_sol_from_rmp(const double*    columns,
                                                const GPtrArray* schedule_sets,
                                                int              num_columns) {
    auto& table = decision_diagram->getDiagram().privateEntity();
    for (auto i = decision_diagram->topLevel(); i >= 0; i--) {
        for (auto& it : table[i]) {
            it.reset_lp_x();
        }
    }

    set_is_integer_solution(true);
    for (int i = 0; i < num_columns; ++i) {
        if (columns[i] > 1e-6) {
            if (columns[i] < 1.0 - 1e-8) {
                set_is_integer_solution(false);
            }

            auto counter = 0u;
            auto tmp = (ScheduleSet*)g_ptr_array_index(schedule_sets, i);
            auto tmp_nodeid(decision_diagram->root());

            while (tmp_nodeid > 1) {
                Job* tmp_j = nullptr;

                if (counter < tmp->job_list->len) {
                    tmp_j = (Job*)g_ptr_array_index(tmp->job_list, counter);
                }

                auto& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    tmp_node.lp_x[1] += columns[i];
                    tmp_nodeid = tmp_node.branch[1];
                    counter++;
                } else {
                    tmp_node.lp_x[0] += columns[i];
                    tmp_nodeid = tmp_node.branch[0];
                }
            }

            assert(tmp_nodeid == 1);
        }
    }

    lp_sol.clear();
    for (auto i = decision_diagram->topLevel(); i > 0; i--) {
        for (auto& it : table[i]) {
            auto value = it.lp_x[1];
            if (value > 1e-6) {
                lp_sol.push_back(
                    BddCoeff(it.get_nb_job(), it.get_weight(), 0.0, value));
            }
            value = it.lp_x[0];
            if (value > 1e-6) {
                lp_sol.push_back(BddCoeff(it.get_nb_job(), it.get_weight(), 0.0,
                                          value, -1, false));
            }
        }
    }

    if (get_is_integer_solution()) {
        added_cuts = false;
        std::cout << "FOUND INTEGER SOLUTION"
                  << "\n";
        getchar();
    } else {
        auto generator =
            ZeroHalfCuts(convex_constr_id, convex_rhs, &reformulation_model,
                         decision_diagram->root(), &table);

        generator.generate_cuts();

        added_cuts = false;
        auto cut_list = generator.get_cut_list();
        for (auto& it : cut_list) {
            reformulation_model.add_constraint(std::move(it));
        }

        added_cuts = (cut_list.size() > 0);
    }

    ColorWriterEdgeX  edge_writer(mip_graph, &table);
    ColorWriterVertex vertex_writer(mip_graph, table);
    auto              file_name = "lp_solution_" + problem_name + "_" +
                     std::to_string(convex_rhs) + ".gv";
    std::ofstream outf(file_name);
    boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    outf.close();
}

void PricerSolverBdd::project_solution(Solution* sol) {
    auto& table = decision_diagram->getDiagram().privateEntity();
    // double*            x = new double[num_edges(mip_graph)]{};
    std::fill(solution_x.get(), solution_x.get() + get_nb_edges(), 0.0);

    for (int i = 0; i < sol->nb_machines; ++i) {
        auto counter = 0u;
        auto tmp = sol->part[i].machine;
        auto tmp_nodeid(decision_diagram->root());

        while (tmp_nodeid > 1) {
            Job* tmp_j = nullptr;

            if (counter < tmp->len) {
                tmp_j = (Job*)g_ptr_array_index(tmp, counter);
            }

            auto& tmp_node = table.node(tmp_nodeid);

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
    // auto& table = decision_diagram->getDiagram().privateEntity();
    // ColorWriterEdgeX   edge_writer(mip_graph, solution_x.get());
    // ColorWriterVertex  vertex_writer(mip_graph, table);
    // string             file_name =
    //     "solution_" + problem_name + "_" + std::to_string(num_machines) +
    //     ".gv";
    // std::ofstream outf(file_name);
    // boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    // outf.close();
}

bool PricerSolverBdd::check_schedule_set(GPtrArray* set) {
    // guint              weight = 0;
    auto& table = decision_diagram->getDiagram().privateEntity();
    auto  tmp_nodeid(decision_diagram->root());
    auto  counter = 0u;

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        Job*  tmp_j = nullptr;

        if (counter < set->len) {
            tmp_j = (Job*)g_ptr_array_index(set, counter);
        }

        if (tmp_j == tmp_node.get_job()) {
            tmp_nodeid = tmp_node.branch[1];
            counter++;
        } else {
            tmp_nodeid = tmp_node.branch[0];
        }
    }

    return (tmp_nodeid == 1 && counter == set->len);
}

void PricerSolverBdd::make_schedule_set_feasible([
    [maybe_unused]] GPtrArray* set) {}

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
