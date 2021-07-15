#include "PricerSolverBdd.hpp"
#include <fmt/core.h>                               // for print
#include <gurobi_c++.h>                             // for GRBLinExpr
#include <algorithm>                                // for find, min, all_of
#include <array>                                    // for array, array<>...
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/graph/adjacency_list.hpp>           // for source, target
#include <boost/graph/detail/adjacency_list.hpp>    // for undirected_edg...
#include <boost/graph/detail/edge.hpp>              // for edge_desc_impl
#include <boost/graph/graphviz.hpp>                 // for write_graphviz
#include <boost/iterator/iterator_facade.hpp>       // for operator!=
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/pending/property.hpp>                  // for lookup_one_pro...
#include <boost/range/irange.hpp>                      // for integer_iterator
#include <cassert>                                     // for assert
#include <cmath>                                       // for fabs
#include <cstddef>                                     // for size_t
#include <ext/alloc_traits.h>                          // for __alloc_traits...
#include <limits>                                      // for numeric_limits
#include <list>                                        // for operator==, list
#include <memory>                                      // for allocator, mak...
#include <ostream>                                     // for operator<<
#include <range/v3/action/action.hpp>                  // for operator|=
#include <range/v3/action/remove_if.hpp>               // for remove_if, rem...
#include <range/v3/algorithm/all_of.hpp>               // for all_of
#include <range/v3/algorithm/fill.hpp>                 // for fill, fill_fn
#include <range/v3/algorithm/max.hpp>                  // for max, max_fn
#include <range/v3/functional/comparisons.hpp>         // for less
#include <range/v3/functional/identity.hpp>            // for identity
#include <range/v3/iterator/diffmax_t.hpp>             // for operator<=
#include <range/v3/iterator/reverse_iterator.hpp>      // for reverse_cursor
#include <range/v3/iterator/unreachable_sentinel.hpp>  // for operator==
#include <range/v3/numeric/partial_sum.hpp>            // for partial_sum
#include <range/v3/view/all.hpp>                       // for all_t
#include <range/v3/view/drop.hpp>                      // for drop, drop_fn
#include <range/v3/view/enumerate.hpp>                 // for enumerate_fn
#include <range/v3/view/iota.hpp>                      // for iota_view, ints
#include <range/v3/view/join.hpp>                      // for join_view, joi...
#include <range/v3/view/partial_sum.hpp>               // for partial_sum
#include <range/v3/view/reverse.hpp>                   // for reverse_view
#include <range/v3/view/subrange.hpp>                  // for subrange
#include <range/v3/view/take.hpp>                      // for take_view, take
#include <range/v3/view/zip.hpp>                       // for zip_view, zip
#include <range/v3/view/zip_with.hpp>                  // for iter_zip_with_...
#include <set>                                         // for operator==
#include <span>                                        // for span
#include <string>                                      // for char_traits
#include <tuple>                                       // for tuple, get
#include <utility>                                     // for pair, move
#include <vector>                                      // for vector, _Bit_r...
#include "CardinalityPaths.hpp"                        // for CardinalityPaths
#include "Column.h"                                    // for ScheduleSet
#include "Instance.h"                                  // for Instance
#include "Job.h"                                       // for Job, value_dif...
#include "Label.hpp"                                   // for Label
#include "MipGraph.hpp"                                // for MipGraph, Colo...
#include "ModelInterface.hpp"                          // for BddCoeff, Refo...
#include "NodeBdd.hpp"                                 // for NodeBdd
#include "NodeBddStructure.hpp"                        // for DdStructure
#include "NodeBddTable.hpp"                            // for NodeTableEntity
#include "NodeData.h"                                  // for NodeData
#include "NodeId.hpp"                                  // for NodeId
#include "OptimalSolution.hpp"                         // for OptimalSolution
#include "PricerConstruct.hpp"                         // for PricerConstruct
#include "PricerSolverBase.hpp"                        // for PricerSolverBa...
#include "Solution.hpp"                                // for Machine, VecJo...
#include "ZeroHalfCuts.hpp"                            // for ZeroHalfCuts
#include "gurobi_c.h"                                  // for GRB_EQUAL, GRB...
#include "lp.h"                                        // for lp_interface_g...
#include "util.h"                                      // for dbg_lvl
#include "util/MyList.hpp"                             // for MyList

PricerSolverBdd::PricerSolverBdd(const Instance& instance)
    : PricerSolverBase(instance),
      decision_diagram(PricerConstruct(instance)),
      size_graph{decision_diagram.size()},
      ordered_jobs_new(instance.vector_pair),
      original_model(reformulation_model),
      H_min{instance.H_min},
      H_max(instance.H_max) {
    // remove_layers_init();

    auto& table = *(decision_diagram.getDiagram());

    auto i = decision_diagram.root().row();
    ordered_jobs_new |=
        ranges::actions::remove_if([&]([[maybe_unused]] const auto& tmp) {
            bool remove = ranges::all_of(
                table[i], [&](const auto& n) { return n[1] == 0; });
            --i;
            return remove;
        });

    if (dbg_lvl() > 0) {
        fmt::print("{0: <{2}}{1}\n", "The new number of layers",
                   ordered_jobs_new.size(), ALIGN);
    }

    decision_diagram.compressBdd();
    init_table();
    cleanup_arcs();
    // check_infeasible_arcs();
    bottum_up_filtering();
    topdown_filtering();
    construct_mipgraph();
    init_coeff_constraints();
    if (dbg_lvl() > 0) {
        fmt::print("ENDING CONSTRUCTION\n\n");
    }
}

PricerSolverBdd::PricerSolverBdd(const PricerSolverBdd& src)
    : PricerSolverBase(src),
      decision_diagram(src.decision_diagram),
      size_graph(src.size_graph),
      nb_removed_edges(src.nb_removed_edges),
      nb_removed_nodes(src.nb_removed_nodes),
      ordered_jobs_new(src.ordered_jobs_new),
      mip_graph(src.mip_graph),
      original_model(src.original_model),
      H_min(src.H_min),
      H_max(src.H_max) {
    // remove_layers_init();
    // decision_diagram.compressBdd();
    // size_graph = decision_diagram.size();
    // init_table();
    cleanup_arcs();
    // check_infeasible_arcs();
    bottum_up_filtering();
    topdown_filtering();
    construct_mipgraph();
    // init_coeff_constraints();
}

PricerSolverBdd::~PricerSolverBdd() = default;

void PricerSolverBdd::construct_mipgraph() {
    mip_graph.clear();
    auto& table = *(decision_diagram.getDiagram());

    auto index{0U};
    for (auto&& [i, row] :
         table | ranges::views::take(decision_diagram.root().row() + 1) |
             ranges::views::enumerate | ranges::views::reverse) {
        for (auto&& [j, it] : row | ranges::views::enumerate) {
            if (NodeId(i, j) != 0 && (it.calc[1] || it.calc[0])) {
                it.key = boost::add_vertex({index++, NodeId(i, j)}, mip_graph);
            }
        }
    }

    auto edge_index{0U};

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (it[0] != 0 && it.calc[0]) {
            auto& n0 = table.node(it[0]);
            add_edge(it.key, n0.key, {edge_index++, false, GRBVar()},
                     mip_graph);
        }

        if (it[1] != 0 && it.calc[1]) {
            auto& n1 = table.node(it[1]);
            add_edge(it.key, n1.key, {edge_index++, true, GRBVar()}, mip_graph);
        }
    }
}

void PricerSolverBdd::init_coeff_constraints() {
    auto& table = *(decision_diagram.getDiagram());
    original_model.clear_all_coeff();

    for (auto i = decision_diagram.root().row(); i > 0UL; i--) {
        for (auto& it : table[i]) {
            it.add_coeff_list_clear();
            for (const auto&& [c, constr] :
                 original_model | ranges::views::take(convex_constr_id) |
                     ranges::views::enumerate) {
                VariableKeyBase key_aux(it.get_nb_job(), it.get_weight());
                auto            coeff = (*constr.get_constr())(key_aux);
                if (fabs(coeff) > EPS_SOLVER) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff, 0.0, c)};
                    constr.add_coeff_to_list(ptr_coeff);
                    it.add_coeff_list(ptr_coeff, true);
                }
            }
        }
    }

    auto&           root = decision_diagram.root();
    auto&           root_node = table.node(root);
    VariableKeyBase key_aux(root_node.get_nb_job(), root_node.get_weight(),
                            true);
    auto*           constr = reformulation_model[convex_constr_id].get();
    auto            coeff = (*constr)(key_aux);
    if (fabs(coeff) > EPS_SOLVER) {
        auto ptr_coeff_high{std::make_shared<BddCoeff>(
            root_node.get_nb_job(), root_node.get_weight(), coeff, true, true)};
        original_model[convex_constr_id].add_coeff_to_list(ptr_coeff_high);
        auto ptr_coeff_low{std::make_shared<BddCoeff>(root_node.get_nb_job(),
                                                      root_node.get_weight(),
                                                      coeff, false, true)};
        original_model[convex_constr_id].add_coeff_to_list(ptr_coeff_low);
    }
}

void PricerSolverBdd::update_coeff_constraints() {
    auto nb_constr = original_model.get_nb_constraints();

    for (auto& j : reformulation_model | ranges::views::drop(nb_constr)) {
        original_model.add_constraint(j);
    }

    auto& table = *(decision_diagram.getDiagram());
    for (auto i = decision_diagram.root().row(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (auto&& [c, constr] : original_model |
                                          ranges::views::drop(nb_constr) |
                                          ranges::views::enumerate) {
                VariableKeyBase key_high{it.get_nb_job(), it.get_weight(),
                                         true};

                auto coeff_high = (*constr.get_constr())(key_high);
                if (fabs(coeff_high) > EPS_SOLVER) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff_high, 0.0,
                        nb_constr + c)};
                    constr.add_coeff_to_list(ptr_coeff);
                    it.add_coeff_list(ptr_coeff, true);
                }

                BddCoeff key_low(it.get_nb_job(), it.get_weight(), 0.0, 0.0,
                                 nb_constr + c, false);
                auto     coeff_low = (*constr.get_constr())(key_low);
                if (fabs(coeff_low) > EPS_SOLVER) {
                    std::shared_ptr<BddCoeff> ptr_coeff{
                        std::make_shared<BddCoeff>(it.get_nb_job(),
                                                   it.get_weight(), coeff_low,
                                                   0.0, nb_constr + c, false)};
                    constr.add_coeff_to_list(ptr_coeff);
                    it.add_coeff_list(ptr_coeff, false);
                }
            }
        }
    }
}

void PricerSolverBdd::init_table() {
    auto& table = *(decision_diagram.getDiagram());
    /** init table */
    auto& root = table.node(decision_diagram.root());
    root.init_node(0, true);
    root.set_node_id_label(decision_diagram.root());
    root.all = boost::dynamic_bitset<>{convex_constr_id, 0};

    for (auto i : ranges::views::ints(0UL, decision_diagram.root().row() + 1) |
                      ranges::views::reverse) {
        for (auto it = 0UL; it < table[i].size(); it++) {
            if (i != 0UL) {
                auto  layer = ordered_jobs_new.size() - i;
                auto& tmp_pair = ordered_jobs_new[layer];
                auto  node_id = NodeId(i, it);
                auto& node = table.node(node_id);
                auto* aux_job = tmp_pair.first;
                auto  w = node.get_weight();
                auto  p = aux_job->processing_time;

                auto& n0 = table.node(node[0]);
                auto& n1 = table.node(node[1]);
                n0.set_node_id_label(node[0]);
                n1.set_node_id_label(node[1]);
                node.set_job_label(aux_job);

                node.ptr_node_id = NodeId(i, it);
                node.set_job(aux_job);
                n0.init_node(w);
                n1.init_node(w + p);
                node.cost = {0.0, static_cast<double>(
                                      aux_job->weighted_tardiness_start(w))};

                n0.in_degree[0]++;
                // n0.in_edges[0].push_back(node.ptr_node_id);
                n1.in_degree[1]++;
                // n1.in_edges[1].push_back(node.ptr_node_id);

            } else {
                auto& node = table.node(NodeId(i, it));
                node.set_job(nullptr);
            }
        }
    }
}

void PricerSolverBdd::insert_constraints_lp(NodeData* pd) {
    lp_interface_get_nb_rows(pd->RMP.get(), &(pd->nb_rows));
    auto nb_new_constraints = reformulation_model.size() - pd->nb_rows;

    fmt::print("nb rows initial {} {} {}\n", pd->nb_rows,
               reformulation_model.size(), nb_new_constraints);

    assert((nb_new_constraints <=
            (pd->id_pseudo_schedules - pd->id_next_var_cuts)));
    std::vector<int>    starts(nb_new_constraints + 1);
    std::vector<char>   sense(nb_new_constraints);
    std::vector<double> rhs(nb_new_constraints);
    std::vector<int>    column_ind;
    std::vector<double> coeff;

    int pos = 0;
    for (auto&& [c, constr] : reformulation_model |
                                  ranges::views::drop(pd->nb_rows) |
                                  ranges::views::enumerate) {
        // auto* constr = reformulation_model[pd->nb_rows + c].get();

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

        for (auto&& [i, it] : pd->localColPool | ranges::views::enumerate) {
            auto& table = *(decision_diagram.getDiagram());
            auto  tmp_nodeid(decision_diagram.root());

            auto coeff_val = 0.0;
            auto job_it = it->job_list.begin();
            while (tmp_nodeid > 1) {
                auto& tmp_node = table.node(tmp_nodeid);
                auto* tmp_j =
                    (job_it != it->job_list.end()) ? *job_it : nullptr;

                VariableKeyBase key(tmp_node.get_nb_job(),
                                    tmp_node.get_weight(),
                                    tmp_j == tmp_node.get_job());
                coeff_val += (*constr)(key);
                if (key.get_high()) {
                    tmp_nodeid = tmp_node[1];
                    ++job_it;
                } else {
                    tmp_nodeid = tmp_node[0];
                }
            }

            assert((tmp_nodeid == 1));

            if (fabs(coeff_val) > EPS_SOLVER) {
                column_ind.push_back(pd->id_pseudo_schedules +
                                     static_cast<int>(i));
                coeff.push_back(coeff_val);
                pos++;
            }
        }
    }

    starts[nb_new_constraints] = pos;

    lp_interface_addrows(pd->RMP.get(), static_cast<int>(nb_new_constraints),
                         static_cast<int>(coeff.size()), starts.data(),
                         column_ind.data(), coeff.data(), sense.data(),
                         rhs.data(), nullptr);
    lp_interface_get_nb_rows(pd->RMP.get(), &(pd->nb_rows));
    pd->pi.resize(reformulation_model.size(), 0.0);
    pd->slack.resize(reformulation_model.size(), 0.0);
    pd->rhs.resize(reformulation_model.size(), 0.0);
    lp_interface_get_rhs(pd->RMP.get(), (pd->rhs).data());
    pd->lhs_coeff.resize(reformulation_model.size(), 0.0);
    pd->id_row.resize(reformulation_model.size(), 0.0);
    pd->coeff_row.resize(reformulation_model.size(), 0.0);
}

double PricerSolverBdd::compute_reduced_cost(const OptimalSolution<>& sol,
                                             double*                  pi,
                                             double*                  lhs) {
    double    result = sol.cost;
    auto&     table = *decision_diagram.getDiagram();
    auto      tmp_nodeid(decision_diagram.root());
    auto      it = sol.jobs.begin();
    std::span aux_lhs{lhs, reformulation_model.size()};
    std::span aux_pi{pi, reformulation_model.size()};
    ranges::fill(aux_lhs, 0.0);

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = (it != sol.jobs.end()) ? *it : nullptr;

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node[1];
            ++it;
            auto* constr = reformulation_model[key.get_j()].get();
            auto  dual = aux_pi[key.get_j()];
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
                aux_lhs[key.get_j()] += coeff;
            }
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             ++c) {
            auto* constr = reformulation_model[c].get();
            auto  dual = aux_pi[c];
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
                aux_lhs[c] += coeff;
            }
        }
    }

    auto*           constr = reformulation_model[convex_constr_id].get();
    auto            dual = aux_pi[convex_constr_id];
    VariableKeyBase k(0, 0, true);
    auto            coeff = (*constr)(k);
    result -= coeff * dual;
    aux_lhs[convex_constr_id] += coeff;

    return result;
}

double PricerSolverBdd::compute_subgradient(const OptimalSolution<>& sol,
                                            double* sub_gradient) {
    double    result = sol.cost;
    auto&     table = *decision_diagram.getDiagram();
    auto      tmp_nodeid(decision_diagram.root());
    auto      it = sol.jobs.begin();
    auto      nb_constraints = reformulation_model.size();
    auto      rhs = -reformulation_model[convex_constr_id]->get_rhs();
    std::span aux_subgradient{sub_gradient, nb_constraints};

    for (auto&& [i, constr] : reformulation_model | ranges::views::enumerate) {
        aux_subgradient[i] = constr->get_rhs();
    }

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = it != sol.jobs.end() ? *it : nullptr;

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node[1];
            ++it;
            auto* constr = reformulation_model[key.get_j()].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_subgradient[key.get_j()] -= coeff * rhs;
            }
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            // auto dual = pi[c];
            auto* constr = reformulation_model[c].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_subgradient[c] -= coeff * rhs;
            }
        }
    }

    aux_subgradient[convex_constr_id] += rhs;
    assert(aux_subgradient[convex_constr_id] == 0.0);
    // sub_gradient[nb_jobs] = 0.0;
    assert(tmp_nodeid == 1);

    return result;
}

bool PricerSolverBdd::refinement_structure(
    const std::vector<std::shared_ptr<Column>>& paths) {
    bool  refined_structure = false;
    auto& table = *decision_diagram.getDiagram();
    for (auto& path : paths) {
        std::vector<NodeId> P;
        std::vector<bool>   L;
        std::vector<NodeId> S;
        std::vector<size_t> index_P;

        S.resize(jobs.size(), NodeId());
        index_P.resize(jobs.size(), 0UL);

        P.push_back(decision_diagram.root());
        auto job_it = path->job_list.begin();
        Job* conflict_job{nullptr};

        while (P.back() > 1) {
            auto& tmp_node{table.node(P.back())};
            auto* tmp_job{job_it != path->job_list.end() ? *job_it : nullptr};

            if (tmp_job == tmp_node.get_job()) {
                if (S[tmp_job->job] != 0) {
                    conflict_job = tmp_job;
                    break;
                }
                S[tmp_job->job] = P.back();
                index_P[tmp_job->job] = P.size() - 1;
                P.push_back(tmp_node[1]);
                L.push_back(true);
                ++job_it;
            } else {
                P.push_back(tmp_node[0]);
                L.push_back(false);
            }
        }

        if (conflict_job) {
            auto nodeid_new = P[index_P[conflict_job->job]];
            for (auto&& label :
                 L | ranges::views::drop(index_P[conflict_job->job])) {
                auto& p = table.node(nodeid_new);
                auto& c = table.node(p[label]);
                auto  w = c;
                nodeid_new =
                    NodeId(p[label].row(), table[p[label].row()].size());

                if (conflict_job == c.get_job()) {
                    w[0] = c[0];
                    w[1] = 0;
                }

                w.set_node_id_label(nodeid_new);
                w.all = boost::dynamic_bitset<>{convex_constr_id, 0};

                table[p[label].row()].emplace_back(w);

                p[label] = nodeid_new;
            }
            refined_structure = true;
        }
    }

    if (refined_structure) {
        decision_diagram.compressBdd();
        nb_removed_nodes -= size_graph;
        size_graph = decision_diagram.size();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }

    return refined_structure;
}

void PricerSolverBdd::enumerate_columns() {
    auto& table = *decision_diagram.getDiagram();
    auto  cursor = 0;
    auto  begin = true;
    std::vector<std::tuple<NodeId, bool, boost::dynamic_bitset<>>> path{};
    auto iterations = boost::multiprecision::cpp_int{};
    auto empty = boost::dynamic_bitset<>{jobs.size(), 0};

    do {
        auto f = begin ? decision_diagram.root() : NodeId(0, 0);
        for (;;) {
            begin = false;
            while (f > 1) { /* down */
                auto const& s = table[f.row()][f.col()];
                auto& set = !path.empty() ? std::get<2>(path.back()) : empty;

                if (s[0] != 0) {
                    cursor = path.size();
                    path.emplace_back(f, false, set);
                    f = s[0];
                } else if (!set[s.get_nb_job()]) {
                    path.emplace_back(f, true, set);
                    std::get<2>(path.back())[s.get_nb_job()] = true;
                    f = s[1];
                } else {
                    f = 0;
                }
            }

            if (f == 1) {
                ++iterations;
                break; /* found */
            }

            for (; cursor >= 0; --cursor) { /* up */
                auto&       sel = path[cursor];
                auto const& ss =
                    table[std::get<0>(sel).row()][std::get<0>(sel).col()];
                if (!std::get<1>(sel) && ss[1] != 0 &&
                    !std::get<2>(sel)[ss.get_nb_job()]) {
                    f = std::get<0>(sel);
                    std::get<1>(sel) = true;
                    path.resize(cursor + 1);
                    std::get<2>(path.back()) = std::get<2>(sel);
                    std::get<2>(path.back())[ss.get_nb_job()] = true;
                    f = decision_diagram.child(f, 1);
                    break;
                }
            }

            if (cursor < 0) { /* end() state */
                fmt::print("number of elementary paths {}\n", iterations.str());
                return;
            }
        }
    } while (true);
}

void PricerSolverBdd::enumerate_columns(double* _pi) {
    auto& table = *get_decision_diagram().getDiagram();
    auto  root_id = get_decision_diagram().root();
    compute_labels(_pi);
    auto reduced_cost = table.node(root_id).backward_label[0].get_f();

    auto cursor = 0;
    auto begin = true;

    std::vector<std::tuple<NodeId, bool, boost::dynamic_bitset<>>> path{};
    auto iterations = boost::multiprecision::cpp_int{};
    auto empty = boost::dynamic_bitset<>{jobs.size(), 0};
    auto aux_nb_machines = static_cast<double>(convex_rhs - 1);

    do {
        auto f = begin ? root_id : NodeId(0, 0);
        for (;;) {
            begin = false;
            while (f > 1) { /* down */
                auto& s = table[f.row()][f.col()];
                auto& set = !path.empty() ? std::get<2>(path.back()) : empty;
                auto  rc = (constLB + aux_nb_machines * reduced_cost +
                           evaluate_rc_arc(s)) < UB;

                if (s[0] != 0) {
                    cursor = path.size();
                    path.emplace_back(f, false, set);
                    f = s[0];
                } else if (!set[s.get_nb_job()] && rc) {
                    path.emplace_back(f, true, set);
                    std::get<2>(path.back())[s.get_nb_job()] = true;
                    f = s[1];
                } else {
                    f = 0;
                }
            }

            if (f == 1) {
                ++iterations;
                break; /* found */
            }

            for (; cursor >= 0; --cursor) { /* up */
                auto& sel = path[cursor];
                auto& ss =
                    table[std::get<0>(sel).row()][std::get<0>(sel).col()];
                auto rc = (constLB + aux_nb_machines * reduced_cost +
                           evaluate_rc_arc(ss)) < UB;
                if (!std::get<1>(sel) && ss[1] != 0 &&
                    !std::get<2>(sel)[ss.get_nb_job()] && rc) {
                    f = std::get<0>(sel);
                    std::get<1>(sel) = true;
                    path.resize(cursor + 1);
                    std::get<2>(path.back()) = std::get<2>(sel);
                    std::get<2>(path.back())[ss.get_nb_job()] = true;
                    f = get_decision_diagram().child(f, 1);
                    break;
                }
            }

            if (cursor < 0) { /* end() state */
                fmt::print("number of elementary paths with rc {}\n",
                           iterations.str());
                return;
            }
        }
    } while (true);
}

double PricerSolverBdd::compute_lagrange(const OptimalSolution<>&   sol,
                                         const std::vector<double>& pi) {
    double result = sol.cost;
    auto   dual_bound = 0.0;

    auto& table = *decision_diagram.getDiagram();
    auto  tmp_nodeid(decision_diagram.root());

    auto it = sol.jobs.begin();
    while (tmp_nodeid > 1) {
        auto&           tmp_node = table.node(tmp_nodeid);
        auto*           tmp_j = it != sol.jobs.end() ? *it : nullptr;
        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            auto dual = pi[key.get_j()];
            auto coeff = (*reformulation_model[key.get_j()])(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
            }

            ++it;
            tmp_nodeid = tmp_node[1];
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto&& [constr, pi_tmp] :
             ranges::views::zip(reformulation_model, pi) |
                 ranges::views::drop(convex_constr_id + 1)) {
            auto coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * pi_tmp;
            }
        }
    }

    result = std::min(0.0, result);

    for (const auto&& [c, constr] :
         reformulation_model | ranges::views::enumerate) {
        if (c == convex_constr_id) {
            continue;
        }

        dual_bound += constr->get_rhs() * pi[c];
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

void PricerSolverBdd::remove_layers_init() {
    auto& table = *(decision_diagram.getDiagram());

    auto i = decision_diagram.root().row();
    ordered_jobs_new |=
        ranges::actions::remove_if([&]([[maybe_unused]] const auto& tmp) {
            bool remove = std::ranges::all_of(
                table[i], [&](const auto& n) { return n[1] == 0; });
            --i;
            return remove;
        });

    if (dbg_lvl() > 0) {
        fmt::print("{0: <{2}}{1}\n", "The new number of layers",
                   ordered_jobs_new.size(), ALIGN);
    }
}

void PricerSolverBdd::remove_layers() {
    // auto& table = *(decision_diagram.getDiagram());

    // auto i = decision_diagram.topLevel();

    // ordered_jobs_new |= ranges::actions::remove_if([&](const auto& tmp) {
    //     auto remove = true;

    //     for (auto& iter : table[i]) {
    //         if (iter.calc[1]) {
    //             remove = false;
    //         } else {
    //             auto& cur_node_1 = iter[1];
    //             cur_node_1 = 0;
    //         }
    //     }
    //     --i;
    //     return remove;
    // });

    // if (dbg_lvl() > 0) {
    //     fmt::print("{0: <{2}}{1}\n", "The new number of layers",
    //                ordered_jobs_new.size(), ALIGN);
    // }
}

void PricerSolverBdd::remove_edges() {
    auto& table = *(decision_diagram.getDiagram());

    /** remove the unnecessary nodes of the bdd */
    for (auto& iter : table | ranges::views::drop(1) | ranges::views::join) {
        if (!iter.calc[1]) {
            auto& cur_node_1 = iter[1];
            iter.ptr_node_id = 0;
            cur_node_1 = 0;
        }

        if (!iter.calc[0]) {
            auto& cur_node_0 = iter[0];
            cur_node_0 = 0;
        }
    }

    decision_diagram.compressBdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram.size();
}

[[maybe_unused]] void PricerSolverBdd::print_representation_file() {
    auto& table = *(decision_diagram.getDiagram());
    // auto  vertex_nodeid_list(get(boost::vertex_name_t(), mip_graph));
    // auto  edge_type_list(get(boost::edge_weight_t(), mip_graph));
    // auto  edge_index_list(get(boost::edge_index_t(), mip_graph));
    auto index_edge{
        std::vector<std::vector<int>>(convex_constr_id, std::vector<int>())};

    auto outfile_file_mip_str =
        problem_name + "_" + std::to_string(convex_rhs) + ".txt";
    std::ofstream out_file_mip(outfile_file_mip_str);

    out_file_mip << boost::num_vertices(mip_graph) << " "
                 << boost::num_edges(mip_graph) << " " << convex_constr_id
                 << " " << convex_rhs << "\n\n";

    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto& head =
            table.node(mip_graph[source(*it.first, mip_graph)].node_id);
        auto& n = table.node(mip_graph[target(*it.first, mip_graph)].node_id);
        if (mip_graph[*it.first].high) {
            double cost =
                head.get_job()->weighted_tardiness_start(head.get_weight());
            out_file_mip << head.key << " " << n.key << " " << cost << '\n';
            index_edge[head.get_nb_job()].push_back(mip_graph[*it.first].id);
        } else {
            out_file_mip << head.key << " " << n.key << " " << 0.0 << '\n';
        }
    }

    out_file_mip << '\n';

    for (auto i = 0UL; i < convex_constr_id; i++) {
        out_file_mip << index_edge[i].size() << " ";
        for (auto& it : index_edge[i]) {
            out_file_mip << it << " ";
        }
        out_file_mip << '\n';
    }

    out_file_mip << '\n';

    for (auto it = vertices(mip_graph); it.first != it.second; it.first++) {
        auto const& node_id = mip_graph[*it.first].node_id;
        auto&       n = table.node(node_id);
        if (node_id > 1) {
            out_file_mip << n.get_nb_job() << " " << n.get_weight() << '\n';
        }
    }

    out_file_mip << "99 99\n";

    out_file_mip.close();
}

bool PricerSolverBdd::evaluate_nodes(double* pi) {
    auto& table = *get_decision_diagram().getDiagram();
    compute_labels(pi);
    auto reduced_cost =
        table.node(get_decision_diagram().root()).backward_label[0].get_f();

    auto removed_edges = false;
    auto nb_removed_edges_evaluate = 0;
    auto aux_nb_machines = static_cast<double>(convex_rhs - 1);

    /** check for each node the Lagrangian dual */
    // for (int i = get_decision_diagram().topLevel(); i > 0; i--) {
    //     for (auto& it : table[i]) {
    for (auto& it :
         table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
             ranges::views ::drop(1) | ranges::views::reverse |
             ranges::views::join) {
        if (((constLB + aux_nb_machines * reduced_cost + evaluate_rc_arc(it) >
              UB - 1.0 + RC_FIXING) ||
             (it.get_job()->weighted_tardiness_start(it.get_weight()) >
              UB - 1.0 + RC_FIXING)) &&
            (it.calc[1])) {
            it.calc[1] = false;
            removed_edges = true;
            add_nb_removed_edges();
            nb_removed_edges_evaluate++;
        }
        if (((constLB + aux_nb_machines * reduced_cost +
                  evaluate_rc_low_arc(it) >
              UB - 1.0 + RC_FIXING)) &&
            (it.calc[0])) {
            it.calc[0] = false;
            removed_edges = true;
            add_nb_removed_edges();
            nb_removed_edges_evaluate++;
        }
    }

    if (removed_edges) {
        if (dbg_lvl() > 0) {
            fmt::print("{0: <{2}}{1}\n",
                       "Number of edges removed by evaluate "
                       "nodes",
                       nb_removed_edges_evaluate, ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                       get_nb_removed_edges(), ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Number of edges", get_nb_edges(),
                       ALIGN);
        }
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }

    return removed_edges;
}

void PricerSolverBdd::add_inequality(std::vector<int> v1, std::vector<int> v2) {
    GRBLinExpr expr1;
    // EdgeVarAccessor   edge_var_list(get(boost::edge_weight2_t(),
    // mip_graph)); EdgeIndexAccessor
    // edge_index_list(get(boost::edge_index_t(), mip_graph));
    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto* data = static_cast<EdgeData*>(it.first->get_property());
        if (std::find(v1.begin(), v1.end(), mip_graph[*it.first].id) !=
            v1.end()) {
            expr1 += mip_graph[*it.first].x;
        }

        if (std::find(v2.begin(), v2.end(), data->id) != v2.end()) {
            expr1 -= mip_graph[*it.first].x;
        }
    }
    model.addConstr(expr1, GRB_EQUAL, 0);
}

void PricerSolverBdd::add_inequality(std::vector<int> v1) {
    GRBLinExpr expr1;
    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto* data = static_cast<EdgeData*>(it.first->get_property());
        if (std::find(v1.begin(), v1.end(), data->id) != v1.end()) {
            expr1 += mip_graph[*it.first].x;
        }
    }
    model.addConstr(expr1, GRB_EQUAL, convex_rhs);
}

void PricerSolverBdd::build_mip() {
    // GRBModel m(*env);
    // model.reset(1);
    // model.set(GRB_IntParam_OutputFlag, 1);
    try {
        fmt::print("Building Mip model for the extended formulation:\n");
        auto& table = *(decision_diagram.getDiagram());

        /** Constructing variables */
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            // auto* data =
            // static_cast<EdgeData*>(it.first->get_property());
            auto& high = mip_graph[*it.first].high;
            auto& x = mip_graph[*it.first].x;
            auto& n =
                table.node(mip_graph[source(*it.first, mip_graph)].node_id);
            if (high) {
                double cost =
                    n.get_job()->weighted_tardiness_start(n.get_weight());
                x = model.addVar(0.0, 1.0, cost, GRB_BINARY);
            } else {
                x = model.addVar(0.0, static_cast<double>(convex_rhs), 0.0,
                                 GRB_INTEGER);
            }
        }

        model.update();

        /** Assignment constraints */
        auto assignment{
            std::vector<GRBLinExpr>(convex_constr_id, GRBLinExpr())};
        auto sense{std::vector<char>(convex_constr_id, GRB_EQUAL)};
        auto rhs{std::vector<double>(convex_constr_id, 1.0)};

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto& high = mip_graph[*it.first].high;
            if (high) {
                auto& x = mip_graph[*it.first].x;
                auto& n =
                    table.node(mip_graph[source(*it.first, mip_graph)].node_id);
                assignment[n.get_nb_job()] += x;
            }
        }

        std::unique_ptr<GRBConstr> assignment_constrs(
            model.addConstrs(assignment.data(), sense.data(), rhs.data(),
                             nullptr, convex_constr_id));
        model.update();

        /** Flow constraints */
        auto num_vertices = boost::num_vertices(mip_graph);
        auto flow_conservation_constr{
            std::vector<GRBLinExpr>(num_vertices, GRBLinExpr())};
        auto sense_flow{std::vector<char>(num_vertices, GRB_EQUAL)};
        auto rhs_flow(std::vector<double>(num_vertices, 0));

        for (auto it = vertices(mip_graph); it.first != it.second; ++it.first) {
            const auto node_id = mip_graph[*it.first].node_id;
            const auto vertex_key = mip_graph[*it.first].index;
            sense_flow[vertex_key] = GRB_EQUAL;
            auto out_edges_it = boost::out_edges(*it.first, mip_graph);

            for (; out_edges_it.first != out_edges_it.second;
                 ++out_edges_it.first) {
                flow_conservation_constr[vertex_key] -=
                    mip_graph[*out_edges_it.first].x;
            }

            auto in_edges_it = boost::in_edges(*it.first, mip_graph);

            for (; in_edges_it.first != in_edges_it.second;
                 ++in_edges_it.first) {
                flow_conservation_constr[vertex_key] +=
                    mip_graph[*in_edges_it.first].x;
            }

            if (node_id == decision_diagram.root()) {
                rhs_flow[vertex_key] = -static_cast<double>(convex_rhs);
            } else if (node_id == 1) {
                rhs_flow[vertex_key] = static_cast<double>(convex_rhs);
            } else {
                rhs_flow[vertex_key] = 0.0;
            }
        }

        std::unique_ptr<GRBConstr> flow_constrs(
            model.addConstrs(flow_conservation_constr.data(), sense_flow.data(),
                             rhs_flow.data(), nullptr, num_vertices));
        model.update();
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            // auto* data =
            // static_cast<EdgeData*>(it.first->get_property());
            auto& high = mip_graph[*it.first].high;
            auto& x = mip_graph[*it.first].x;
            auto& n =
                table.node(mip_graph[source(*it.first, mip_graph)].node_id);
            if (high) {
                x.set(GRB_DoubleAttr_Start, n.best_sol_x[1]);

            } else {
                x.set(GRB_DoubleAttr_Start, n.best_sol_x[0]);
            }
        }
        model.optimize();
    } catch (GRBException& e) {
        fmt::print("Error code = {}\n", e.getErrorCode());
        fmt::print("{}", e.getMessage());
    } catch (...) {
        fmt::print("Exception during optimization\n");
    }
}

void PricerSolverBdd::cleanup_arcs() {
    NodeTableEntity<>& table = *(decision_diagram.getDiagram());

    table.node(0).backward_distance[0] = std::numeric_limits<int>::min();
    table.node(0).backward_distance[1] = std::numeric_limits<int>::min();
    table.node(1).backward_distance[0] = 0;
    table.node(1).backward_distance[1] = 0;
    auto removed_edges = false;
    auto nb_edges_removed_tmp = 0;

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::join) {
        // it.calc = {true, true};
        for (size_t i = 0UL; i < 2; ++i) {
            auto& cur_node = table.node(it[i]);
            it.backward_distance.at(i) =
                ranges::max(cur_node.backward_distance) +
                i * it.get_job()->processing_time;
        }
    }
    /** remove the unnecessary nodes of the bdd */
    for (auto& iter : table |
                          ranges::views::take(decision_diagram.topLevel() + 1) |
                          ranges::views::drop(1) | ranges::views::join) {
        // if (iter.get_weight() + iter.backward_distance[0] < H_min &&
        //     iter.branch[0] != 0) {
        //     iter.calc[0] = false;
        //     removed_edges = true;
        //     nb_edges_removed_tmp++;
        //     nb_removed_edges++;
        // }

        if (iter.get_weight() + iter.backward_distance[1] < H_min) {
            iter.calc[1] = false;
            removed_edges = true;
            nb_edges_removed_tmp++;
            nb_removed_edges++;
        }
    }

    if (removed_edges) {
        if (dbg_lvl() > 0) {
            fmt::print("{0: <{2}}{1}\n", "Number of edges removed by clean up",
                       nb_edges_removed_tmp, ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                       get_nb_removed_edges(), ALIGN);
        }
        remove_layers();
        remove_edges();
    }
}

void PricerSolverBdd::topdown_filtering() {
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = *(decision_diagram.getDiagram());
    auto& root = table.node(decision_diagram.root());
    root.init_node(0, true);
    for (auto& it : table |
                        ranges::views::take(decision_diagram.root().row() + 1) |
                        ranges::views::reverse | ranges::views::join) {
        // for (auto& it : table[i]) {
        it.visited = false;
        it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
        // it.calc[1] = true;
        // }
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        auto& n0 = table.node(it[0]);

        if (n0.visited) {
            n0.all &= it.all;
        } else {
            n0.visited = true;
            n0.all = boost::dynamic_bitset<>{convex_constr_id, 0};
            n0.all |= it.all;
        }

        auto& n1 = table.node(it[1]);

        if (n1.visited) {
            if (n1.all[it.get_nb_job()]) {
                n1.all &= it.all;
                n1.all[it.get_nb_job()] = true;
            } else {
                n1.all &= it.all;
            }
        } else {
            n1.all = boost::dynamic_bitset<>{convex_constr_id, 0};
            n1.all |= it.all;
            n1.all[it.get_nb_job()] = true;
            n1.visited = true;
        }
    }

    for (auto i : ranges::views::ints(1UL, decision_diagram.root().row() + 1) |
                      ranges::views::reverse) {
        for (auto& it : table[i]) {
            if (it.all[it.get_nb_job()]) {
                removed_edges = true;
                it.calc[1] = false;
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
        // remove_layers();
        remove_edges();
        cleanup_arcs();
        // init_table();
        // continue;
    }
}

void PricerSolverBdd::bottum_up_filtering() {
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = *(decision_diagram.getDiagram());
    for (auto i : ranges::views::ints(1UL, decision_diagram.root().row() + 1) |
                      ranges::views::reverse) {
        for (auto& it : table[i]) {
            it.visited = false;
            it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
            // it.calc[1] = true;
        }
    }

    table.node(0).all = boost::dynamic_bitset<>{convex_constr_id, 0};
    table.node(1).all = boost::dynamic_bitset<>{convex_constr_id, 0};
    table.node(0).all.flip();

    for (auto i : ranges::views::ints(1UL, decision_diagram.root().row() + 1)) {
        for (auto& it : table[i]) {
            it.all[it.get_nb_job()] = true;
            it.all |= table.node(it[1]).all;
            it.all &= table.node(it[0]).all;
        }
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (table.node(it[1]).all[it.get_nb_job()]) {
            removed_edges = true;
            it.calc[1] = false;
            nb_removed_edges++;
            nb_edges_removed_tmp++;
        }
    }

    if (removed_edges) {
        // std::cout << "removing edges based on bottum-up iteration\n";
        // std::cout << "Number edges removed bottum-up iteration \t\t= "
        //           << nb_edges_removed_tmp << "\n";
        // std::cout << "Number edges removed total \t\t\t= " <<
        // nb_removed_edges
        //           << "\n";
        // remove_layers();
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
    auto& table = *(decision_diagram.getDiagram());

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::reverse | ranges::views::join) {
        it.visited = false;
        it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
        it.calc = {true, true};
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        auto& n0 = table.node(it[0]);
        n0.all |= it.all;
        auto& n1 = table.node(it[1]);
        n1.all[it.get_nb_job()] = true;
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (!it.all.empty() &&
            it.all.find_first() != boost::dynamic_bitset<>::npos) {
            auto index = it.all.find_first();

            auto max = value_diff_Fij(it.get_weight(), it.get_job(),
                                      jobs[index].get());
            // bool index_bool = (index > (size_t)it.get_nb_job());
            while (index != boost::dynamic_bitset<>::npos && max < 0) {
                index = it.all.find_next(index);
                if (index != boost::dynamic_bitset<>::npos) {
                    int a = value_diff_Fij(it.get_weight(), it.get_job(),
                                           jobs[index].get());
                    if (a > max) {
                        max = a;
                    }
                }
            }

            if (max < 0) {
                removed_edges = true;
                it.calc[1] = false;
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("removing edges based on order\n");
        fmt::print("Number edges removed order = {}\n", nb_edges_removed_tmp);
        fmt::print("Number edges removed total = {}\n", nb_removed_edges);
        remove_layers();
        remove_edges();
        cleanup_arcs();
        // init_table();
    }
}

void PricerSolverBdd::equivalent_paths_filtering() {
    /** init table */
    auto removed_edges = false;
    auto nb_edges_removed_tmp = 0;
    // auto  edge_type_list(get(boost::edge_weight_t(), mip_graph));
    auto& table = *(decision_diagram.getDiagram());
    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        it.visited = false;
        it.all = boost::dynamic_bitset<>{convex_constr_id, 0};
        // it.calc = {true, true};
        auto& n0 = table.node(it[0]);
        n0.in_degree[0]++;
        auto& n1 = table.node(it[1]);
        n1.in_degree[1]++;
    }

    std::vector<size_t> vertices;

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (it.in_degree[0] + it.in_degree[0] >= 2) {
            vertices.push_back(it.key);
        }
    }

    for (auto& it : vertices) {
        auto              start_v = mip_graph[it].node_id;
        auto              num_vertices = boost::num_vertices(mip_graph);
        std::list<NodeId> queue;

        queue.push_back(start_v);
        std::vector<bool>                    visited(num_vertices, false);
        std::vector<bool>                    edge_visited(num_vertices, false);
        std::vector<boost::dynamic_bitset<>> all(
            num_vertices, boost::dynamic_bitset(convex_constr_id, 0));
        std::vector<int> C(num_vertices, 0);

        auto& tmp_n = table.node(start_v);
        visited[tmp_n.key];
        auto stop = false;

        while (!queue.empty()) {
            auto currVertex = queue.front();
            queue.pop_front();
            auto iter = boost::in_edges(table.node(currVertex).key, mip_graph);

            for (; iter.first != iter.second; iter.first++) {
                auto adjVertex =
                    mip_graph[source(*iter.first, mip_graph)].node_id;
                auto& n = table.node(adjVertex);
                auto* data = static_cast<EdgeData*>(iter.first->get_property());
                auto  high = data->high;

                if (!visited[n.key]) {
                    visited[n.key] = true;
                    queue.push_back(adjVertex);
                    if (high) {
                        auto& tmp_node = table.node(n[1]);
                        all[n.key] |= all[tmp_node.key];
                        all[n.key][n.get_nb_job()] = true;
                        C[n.key] =
                            C[tmp_node.key] + n.get_job()->weighted_tardiness(
                                                  tmp_node.get_weight());
                        edge_visited[n.key] = true;
                    } else {
                        auto& tmp_node = table.node(n[0]);
                        all[n.key] |= all[tmp_node.key];
                        C[n.key] = C[tmp_node.key];
                    }
                } else {
                    boost::dynamic_bitset<> tmp;
                    int                     tmp_C{};
                    if (high) {
                        auto& tmp_node = table.node(n[1]);
                        tmp = all[tmp_node.key];
                        tmp[n.get_nb_job()] = true;
                        tmp_C =
                            C[tmp_node.key] + n.get_job()->weighted_tardiness(
                                                  tmp_node.get_weight());
                    } else {
                        auto& tmp_node = table.node(n[0]);
                        tmp = all[tmp_node.key];
                        tmp_C = C[tmp_node.key];
                    }

                    if (all[n.key] == tmp) {
                        NodeId cur{};
                        NodeId prev = adjVertex;
                        if (high) {
                            if (tmp_C > C[n.key]) {
                                cur = n[1];
                            } else {
                                cur = n[0];
                            }
                        } else {
                            if (tmp_C > C[n.key]) {
                                cur = n[0];
                            } else {
                                cur = n[1];
                            }
                        }

                        while (cur != start_v) {
                            auto& node = table.node(cur);
                            if (node.in_degree[0] + node.in_degree[1] > 1) {
                                break;
                            }
                            prev = cur;
                            assert(cur != NodeId(0, 0));
                            assert(cur != NodeId(0, 1));
                            if (edge_visited[node.key]) {
                                cur = node[1];
                            } else {
                                cur = node[0];
                            }
                        }

                        auto& node_delete = table.node(prev);
                        if (edge_visited[node_delete.key] && cur == start_v) {
                            node_delete.calc[1] = false;
                            removed_edges = true;
                            nb_edges_removed_tmp++;
                        } else if (cur == start_v) {
                            node_delete.calc[0] = false;
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
        fmt::print(
            "Number of edges removed by equivalent_path_filtering = "
            "{}\nNumber "
            "of edges removed in total = {}\n",
            nb_edges_removed_tmp, nb_removed_edges);

        remove_layers();
        remove_edges();
        cleanup_arcs();
        construct_mipgraph();
    }
}

void PricerSolverBdd::construct_lp_sol_from_rmp(
    const double*                               columns,
    const std::vector<std::shared_ptr<Column>>& schedule_sets) {
    auto& table = *(decision_diagram.getDiagram());
    for (auto& it : table |
                        ranges::views::take(decision_diagram.root().row() + 1) |
                        ranges::views::join) {
        // for (auto& it : table[i]) {
        it.reset_lp_x();
        // }
    }
    std::span<const double> aux_cols{columns, schedule_sets.size()};

    set_is_integer_solution(true);
    for (auto&& [i, x] : aux_cols | ranges::views::enumerate) {
        if (x > EPS_SOLVER) {
            auto* tmp = schedule_sets[i].get();
            auto  it = tmp->job_list.begin();
            auto  tmp_nodeid(decision_diagram.root());

            while (tmp_nodeid > 1) {
                Job*  tmp_j = (it != tmp->job_list.end()) ? *it : nullptr;
                auto& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    tmp_node.lp_x[1] += x;
                    tmp_nodeid = tmp_node[1];
                    ++it;
                } else {
                    tmp_node.lp_x[0] += x;
                    tmp_nodeid = tmp_node[0];
                }
            }

            assert(tmp_nodeid == 1);
        }
    }

    lp_sol.clear();
    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        it.lp_visited = false;
        auto value = it.lp_x[1];
        if (value > EPS_SOLVER) {
            lp_sol.emplace_back(it.get_nb_job(), it.get_weight(), 0.0, value);
            if (value < 1.0 - EPS_SOLVER) {
                set_is_integer_solution(false);
            }
        }

        // value = it.lp_x[0];
        // if (value > EPS_SOLVER) {
        //     lp_sol.emplace_back(it.get_nb_job(), it.get_weight(),
        //     0.0,
        //                         value, -1, false);
        // }
    }

    if (is_integer_solution && dbg_lvl() > 1) {
        fmt::print("FOUND INTEGER SOLUTION\n\n");
    }

    // ColorWriterEdgeX  edge_writer(mip_graph, &table);
    // ColorWriterVertex vertex_writer(mip_graph, table);
    // auto              file_name = "lp_solution_" + problem_name + "_" +
    //                  std::to_string(convex_rhs) + ".gv";
    // std::ofstream outf(file_name);
    // boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    // outf.close();
}

void PricerSolverBdd::project_sol_on_original_variables(const Sol& _sol) {
    auto& table = *(decision_diagram.getDiagram());
    for (auto& m : _sol.machines) {
        auto tmp_nodeid(decision_diagram.root());
        auto it = m.job_list.begin();

        while (tmp_nodeid > 1) {
            auto* tmp_j = (it != m.job_list.end()) ? *it : nullptr;
            auto& tmp_node = table.node(tmp_nodeid);

            if (tmp_j == tmp_node.get_job()) {
                tmp_node.best_sol_x[1] += 1.0;
                tmp_nodeid = tmp_node[1];
                ++it;
            } else {
                tmp_node.best_sol_x[0] += 1.0;
                tmp_nodeid = tmp_node[0];
            }
        }

        assert(tmp_nodeid == 1);
    }
}

std::vector<std::vector<double>>& PricerSolverBdd::calculate_job_time() {
    ranges::fill(x_bar | ranges::views::join, 0.0);
    ranges::fill(z_bar | ranges::views::join, 0.0);

    for (auto& it : lp_sol) {
        if (it.get_high()) {
            x_bar[it.get_j()][it.get_t()] += it.get_value();
        }
    }

    for (auto&& [x_j, z_j] : ranges::views::zip(x_bar, z_bar)) {
        z_j =
            ranges::views::partial_sum(x_j, std::plus<>{}) | ranges::to_vector;
    }

    return x_bar;
}

void PricerSolverBdd::split_job_time(size_t _job, int _time, bool _left) {
    auto& table = *(decision_diagram.getDiagram());
    auto  removed_edges = false;

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (_left) {
            if (it.get_weight() <= _time && it.get_nb_job() == _job) {
                it.calc[1] = false;
                removed_edges = true;
            }
        } else {
            if (it.get_weight() > _time && it.get_nb_job() == _job) {
                it.calc[1] = false;
                removed_edges = true;
            }
        }
    }

    if (removed_edges) {
        decision_diagram.compressBdd();
        remove_layers();
        remove_edges();
        bottum_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();

        if (dbg_lvl() > 0) {
            auto&             table_bis = *(decision_diagram.getDiagram());
            ColorWriterVertex vertex_writer(mip_graph, table_bis);
            auto file_name = "split_solution_" + problem_name + "_" +
                             std::to_string(_job) + "_" +
                             std::to_string(_time) + "_" +
                             std::to_string(_left) + ".gv";
            std::ofstream outf(file_name);
            boost::write_graphviz(outf, mip_graph, vertex_writer);
            outf.close();
        }
    }
}

int PricerSolverBdd::add_constraints() {
    auto& table = *(decision_diagram.getDiagram());
    if (get_is_integer_solution()) {
        bool added_cuts = false;
        return added_cuts;
    } else {
        auto generator =
            ZeroHalfCuts(convex_constr_id, convex_rhs, &reformulation_model,
                         decision_diagram.root(), &table);

        generator.generate_cuts();

        bool added_cuts = false;
        auto cut_list = generator.get_cut_list();
        for (auto& it : cut_list) {
            reformulation_model.emplace_back(std::move(it));
        }

        added_cuts = !cut_list.empty();
        return added_cuts;
    }
}

void PricerSolverBdd::remove_constraints(int first, int nb_del) {
    original_model.delete_constraints(first, nb_del);
    reformulation_model.delete_constraints(first, nb_del);
}

void PricerSolverBdd::update_rows_coeff(size_t first) {
    for (auto k :
         ranges::views::ints(first, original_model.get_nb_constraints())) {
        auto& aux = original_model.get_coeff_list(k);
        for (auto& it : aux) {
            it->set_row(k);
        }
    }
}

bool PricerSolverBdd::check_schedule_set(const std::vector<Job*>& set) {
    auto& table = *(decision_diagram.getDiagram());
    auto  tmp_nodeid(decision_diagram.root());
    auto  it = set.begin();

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = (it != set.end()) ? *it : nullptr;

        if (tmp_j == tmp_node.get_job()) {
            tmp_nodeid = tmp_node[1];
            ++it;
        } else {
            tmp_nodeid = tmp_node[0];
        }
    }

    return (tmp_nodeid == 1 && it == set.end());
}

double PricerSolverBdd::evaluate_rc_low_arc(NodeBdd<>& n) {
    auto& table = *(get_decision_diagram().getDiagram());
    auto& child = table.node(n[0]);
    return n.forward_label[0].get_f() + child.backward_label[0].get_f() +
           n.reduced_cost[0];
}

void PricerSolverBdd::iterate_zdd() {
    DdStructure<NodeBdd<double>>::const_iterator it = decision_diagram.begin();

    for (; it != decision_diagram.end(); ++it) {
        auto i = (*it).begin();

        for (; i != (*it).end(); ++i) {
            fmt::print("{} ", ordered_jobs_new.size() - *i);
        }

        fmt::print("\n");
    }
}

void PricerSolverBdd::create_dot_zdd(const char* name) {
    auto&             table = *(decision_diagram.getDiagram());
    ColorWriterVertex vertex_writer(mip_graph, table);
    std::ofstream     outf(name);
    boost::write_graphviz(outf, mip_graph, vertex_writer);
    outf.close();
}

void PricerSolverBdd::print_number_nodes_edges() {
    fmt::print("removed edges = {}, removed nodes = {}\n", nb_removed_edges,
               nb_removed_nodes);
}

size_t PricerSolverBdd::get_num_remove_nodes() {
    return nb_removed_nodes;
}

size_t PricerSolverBdd::get_num_remove_edges() {
    return nb_removed_edges;
}

size_t PricerSolverBdd::get_nb_edges() {
    return num_edges(mip_graph);
}

size_t PricerSolverBdd::get_nb_vertices() {
    return num_vertices(mip_graph);
}

bool PricerSolverBdd::structure_feasible() {
    return (get_nb_vertices() != 0 &&
            decision_diagram.root().row() >= jobs.size());
}

int PricerSolverBdd::get_num_layers() {
    return decision_diagram.topLevel();
}

boost::multiprecision::cpp_int PricerSolverBdd::print_num_paths() {
    auto evaluator = CardinalityPaths();
    return decision_diagram.evaluate_backward(evaluator);
}
