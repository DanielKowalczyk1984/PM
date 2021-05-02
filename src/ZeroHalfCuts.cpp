#include "ZeroHalfCuts.hpp"
#include <fmt/core.h>
#include <fmt/format.h>
#include <gurobi_c++.h>
#include <algorithm>
#include <memory>
#include <range/v3/all.hpp>
#include <vector>
#include "ModelInterface.hpp"
#include "NodeBddTable.hpp"
#include "NodeId.hpp"

ZeroHalfCuts::ZeroHalfCuts(size_t              _nb_jobs,
                           size_t              _nb_machines,
                           ReformulationModel* _rmp_model,
                           NodeId const&       _root,
                           NodeTableEntity<>*  _table)
    : env(std::make_unique<GRBEnv>()),
      model(std::make_unique<GRBModel>(*env)),
      nb_jobs(_nb_jobs),
      nb_machines(_nb_machines),
      rmp_model(_rmp_model),
      root(_root),
      table(_table),
      jobs_var(nb_jobs) {
    // Limit how many solutions to collect
    // model->set(GRB_IntParam_PoolSolutions, 124);

    // Limit the search space by setting a gap for the worst possible solution
    // that will be accepted
    // model->set(GRB_DoubleParam_PoolGap, 0.10);

    // do a systematic search for the k-best solutions
    model->set(GRB_IntParam_PoolSearchMode, GRB_PRESOLVE_AGGRESSIVE);
    model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    model->set(GRB_DoubleParam_Heuristics, 1 / HALF);
    model->set(GRB_IntParam_MIPFocus, 1);
    model->set(GRB_DoubleParam_TimeLimit, TIMELIMIT);

    generate_model();
}

bool ZeroHalfCuts::add_cuts() {
    return true;
}

void ZeroHalfCuts::generate_model() {
    try {
        for (auto i : ranges::views::ints(0UL, nb_jobs)) {
            jobs_var[i] =
                model->addVar(0.0, 1.0, 0.0, 'B', fmt::format("jobs_{}", i));
        }

        auto& root_node = table->node(root);
        auto& terminal_node = table->node(1);
        terminal_node.lp_visited = false;
        root_node.sigma =
            model->addVar(0.0, 1.0, 0.0, 'B', fmt::format("sigma_root"));
        node_ids.push_back(root);

        dfs(root);

        GRBLinExpr          expr = 0;
        std::vector<double> coeffs(nb_jobs, 1.0);
        auto                m = nb_machines % 2 == 0 ? 0.0 : 1.0;
        q = model->addVar(0.0, GRB_INFINITY, 0.0, 'I');

        expr.addTerms(coeffs.data(), jobs_var.data(), coeffs.size());
        expr += m * root_node.sigma - m * terminal_node.sigma - HALF * q;
        model->addConstr(expr, '=', 1.0);
        model->update();
        init_table();

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "",
                   ALIGN, e.getMessage());
    }
}

void ZeroHalfCuts::init_table() {
    for (auto& it : (*table) | ranges::views::join) {
        it.in_degree = {};
        // it.in_degree_1 = 0;
        it.in_edges[0].clear();
        it.in_edges[1].clear();
        it.coeff_cut = {};
    }

    for (auto& it : *table | ranges::views::take(root.row() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        auto& n0 = table->node(it[0]);
        auto& n1 = table->node(it[1]);

        n0.in_edges[0].push_back(it.ptr_node_id);
        n0.in_degree[0]++;
        n1.in_edges[1].push_back(it.ptr_node_id);
        n1.in_degree[1]++;
    }
}

void ZeroHalfCuts::init_coeff_cut() {
    for (auto& it : node_ids) {
        auto& node = table->node(it);
        init_coeff_node(&node);
    }

    for (auto& iter : node_ids_lift) {
        auto& node = table->node(iter);
        init_coeff_node(&node);
    }
    node_ids_lift.clear();
}

void ZeroHalfCuts::init_coeff_node(NodeBdd<>* node) {
    for (auto k : ranges::views::ints(0UL, 2UL)) {
        node->coeff_cut.at(k) = 0.0;
        for (auto& it : node->in_edges.at(k)) {
            auto aux = it.lock();
            if (aux) {
                auto& aux_node = table->node(*aux);
                aux_node.coeff_cut.at(k) = 0.0;
            }
        }
    }
}

void ZeroHalfCuts::construct_cut() {
    std::unique_ptr<GenericData> data = std::make_unique<GenericData>();

    auto add_coeff_constr = [&](const auto& it) {
        auto& node = table->node(it);

        for (auto k = 0UL; k < 2; k++) {
            auto coeff = floor(node.coeff_cut.at(k) / HALF);
            if (coeff > EPS_CUT) {
                data->add_coeff_hash_table(node.get_nb_job(), node.get_weight(),
                                           k, -coeff);
            }

            for (auto& iter : node.in_edges.at(k)) {
                auto aux = iter.lock();
                if (aux) {
                    auto& aux_node = table->node(*aux);
                    auto  coeff_in = floor(aux_node.coeff_cut.at(k) / 2);
                    if (coeff_in < -EPS_CUT) {
                        data->add_coeff_hash_table(aux_node.get_nb_job(),
                                                   aux_node.get_weight(), k,
                                                   -coeff_in);
                    }
                }
            }
        }
    };

    // auto print_node_ids = [&](const auto& it) {
    //     auto& node = table->node(it);
    //     if (it > 1 && node.sigma.get(GRB_DoubleAttr_Xn) > 0.0) {
    //         fmt::print("node {} {} {} {} {} | ", node.get_nb_job(),
    //                    node.get_weight(), node.sigma.get(GRB_DoubleAttr_Xn),
    //                    node.lp_x[0], node.lp_x[1]);
    //     } else if (node.sigma.get(GRB_DoubleAttr_Xn) > 0.0) {
    //         fmt::print("Terminal Node | ");
    //     }
    // };

    std::ranges::for_each(node_ids, add_coeff_constr);
    // std::for_each(node_ids_lift.begin(), node_ids_lift.end(),
    // add_coeff_constr); std::for_each(node_ids.begin(), node_ids.end(),
    // print_node_ids); fmt::print("\n");

    auto& root_node = table->node(root);
    auto& terminal_node = table->node(1);
    auto  rhs = 0.0;
    rhs += (root_node.sigma.get(GRB_DoubleAttr_Xn) > EPS_CUT)
               ? static_cast<double>(nb_machines)
               : 0.0;
    rhs += (terminal_node.sigma.get(GRB_DoubleAttr_Xn) > EPS_CUT)
               ? -static_cast<double>(nb_machines)
               : 0.0;

    for (auto& it : jobs_var) {
        auto x = it.get(GRB_DoubleAttr_Xn);
        if (x > EPS_CUT) {
            rhs += 1.0;
        }
    }
    // fmt::print("test rhs = {}\n", rhs);
    std::shared_ptr<ConstraintGeneric> constr{
        std::make_shared<ConstraintGeneric>(data.release(),
                                            -floor(rhs / HALF))};
    // data->list_coeff();
    // fmt::print("RHS = {}\n",
    //            -floor((2.0 * q.get(GRB_DoubleAttr_Xn) + 1 + rhs) / 2.0));
    for (auto& it : cut_list) {
        if (*constr == *it) {
            return;
        }
    }
    cut_list.push_back(std::move(constr));
}

std::vector<std::shared_ptr<ConstraintGeneric>> ZeroHalfCuts::get_cut_list() {
    return cut_list;
}

void ZeroHalfCuts::generate_cuts() {
    try {
        model->optimize();
        auto status = model->get(GRB_IntAttr_Status);

        if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE ||
            status == GRB_UNBOUNDED) {
            fmt::print(
                "The model cannot be solved because it is infeasible or "
                "unbounded\n");
            return;
        }

        if (status != GRB_OPTIMAL) {
            fmt::print("Optimization was stopped with status {}\n", status);
        }

        // Print number of solutions stored
        auto nb_solutions = model->get(GRB_IntAttr_SolCount);
        fmt::print("Number of solutions found: {}\n", nb_solutions);

        for (auto i = 0; i < nb_solutions; i++) {
            model->set(GRB_IntParam_SolutionNumber, i);
            init_coeff_cut();

            bool add = false;
            auto calc_coeff_cut = [&](const auto& it) {
                auto& node = table->node(it);
                auto  x = node.sigma.get(GRB_DoubleAttr_Xn);
                if (x > EPS_CUT) {
                    if (it > 1) {
                        node.coeff_cut[0] += 1.0;
                        node.coeff_cut[1] += 1.0;
                        // fmt::print("Node {} {}\n", node.get_nb_job(),
                        //    node.get_weight());
                        add = true;
                    }

                    for (auto j : ranges::views::ints(0UL, 2UL)) {
                        for (auto& it_aux : node.in_edges.at(j)) {
                            auto aux = it_aux.lock();
                            if (aux) {
                                auto& aux_node = table->node(*aux);
                                aux_node.coeff_cut.at(j) -= 1.0;
                            }
                        }
                    }
                }
            };

            for (auto k = root.row(); k > 0; k--) {
                for (auto& it : (*table)[k]) {
                    auto j = it.get_nb_job();
                    auto x = jobs_var[j].get(GRB_DoubleAttr_Xn);
                    if (x > EPS_CUT) {
                        it.coeff_cut[1] += 1.0;
                    }
                }
            }

            // for (auto i = 0; i < nb_jobs; i++) {
            //     auto x = jobs_var[i].get(GRB_DoubleAttr_Xn);
            //     if (x > 1e-4) {
            //         fmt::print("{} ", i);
            //     }
            // }
            // fmt::print("\n");

            std::ranges::for_each(node_ids, calc_coeff_cut);

            // for (auto& iter : node_ids) {
            //     auto& node = table->node(iter);
            //     if (node.branch[0] <= 1) {
            //         continue;
            //     }
            //     auto x = node.sigma.get(GRB_DoubleAttr_Xn);
            //     if ((x > 1e-4 && node.lp_x[0] < 1e-6)) {
            //         auto& child = table->node(node.branch[0]);
            //         for (int j = 0; j < 2; j++) {
            //             child.coeff_cut[j] += 1.0;
            //             // if (j) {
            //             //     auto y = jobs_var[child.get_nb_job()].get(
            //             //         GRB_DoubleAttr_Xn);
            //             //     if (y > 1e-4) {
            //             //         child.coeff_cut[j] += 1.0;
            //             //     }
            //             // }
            //             for (auto& it : child.in_edges[j]) {
            //                 auto aux = it.lock();
            //                 if (aux) {
            //                     auto& aux_node = table->node(*aux);
            //                     aux_node.coeff_cut[j] -= 1.0;
            //                     // auto y =
            //                     jobs_var[aux_node.get_nb_job()].get(
            //                     //     GRB_DoubleAttr_Xn);
            //                     // if (y == 1.0 && !aux_node.lp_visited) {
            //                     //     aux_node.coeff_cut[j] += 1.0;
            //                     // }
            //                 }
            //             }
            //         }
            //         node_ids_lift.push_back(node.branch[0]);
            //         dfs_lift(node.branch[0]);
            //     }
            // }
            if (add) {
                construct_cut();
            }
            // fmt::print("---------------------------------------------\n");
        }
        // getchar();

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "",
                   ALIGN, e.getMessage());
    }
}

void ZeroHalfCuts::dfs(const NodeId& v) {
    auto& node = table->node(v);
    node.lp_visited = true;

    for (auto i : ranges::views::ints(0UL, 2UL)) {
        if (node.lp_x.at(i) > EPS_CUT) {
            auto& child = table->node(node.at(i));
            if (!child.lp_visited) {
                if (node.at(i) == 1) {
                    child.sigma = model->addVar(0.0, 1.0, 0.0, 'B',
                                                fmt::format("sigma_terminal"));
                } else {
                    child.sigma = model->addVar(
                        0.0, 1.0, 0.0, 'B',
                        fmt::format("sigma_{}_{}", child.get_nb_job(),
                                    child.get_weight()));
                }

                node_ids.push_back(node.at(i));
                auto& s_source = node.sigma;
                auto& s_head = child.sigma;
                auto  str_y = fmt::format("y_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto  str_r = fmt::format("r_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto& y = node.y.at(i) = model->addVar(
                    0.0, GRB_INFINITY, node.lp_x.at(i), 'B', str_y);
                auto& r = node.r.at(i) =
                    model->addVar(0.0, GRB_INFINITY, 0.0, 'I', str_r);
                GRBLinExpr expr = -s_head + s_source - y - HALF * r;
                if (i) {
                    expr += jobs_var[node.get_nb_job()];
                }
                model->addConstr(expr, '=', 0.0);
                dfs(node.at(i));
            } else {
                auto& s_source = node.sigma;
                auto& s_head = child.sigma;
                auto  str_y = fmt::format("y_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto  str_r = fmt::format("r_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto& y = node.y.at(i) = model->addVar(
                    0.0, GRB_INFINITY, (node.lp_x).at(i), 'B', str_y);
                auto& r = node.r.at(i) =
                    model->addVar(0.0, GRB_INFINITY, 0.0, 'I', str_r);
                GRBLinExpr expr = -s_head + s_source - y - HALF * r;
                if (i) {
                    expr += jobs_var[node.get_nb_job()];
                }

                model->addConstr(expr, '=', 0.0);
            }
        }
    }
}

void ZeroHalfCuts::dfs_lift(const NodeId& v) {
    auto& node = table->node(v);

    for (auto k : ranges::views::ints(0UL, 2UL)) {
        if (node.at(k) <= 1) {
            continue;
        }
        auto& child_node = table->node(node.at(k));

        if (node.lp_x.at(k) == 0.0 && !child_node.lp_visited) {
            child_node.coeff_cut[0] += 1.0;
            child_node.coeff_cut[1] += 1.0;

            for (auto j : ranges::views::ints(0UL, 2UL)) {
                for (auto& it : child_node.in_edges.at(j)) {
                    auto aux = it.lock();
                    if (aux) {
                        auto& aux_node = table->node(*aux);
                        aux_node.coeff_cut.at(j) -= 1.0;
                    }
                }
            }
            node_ids_lift.push_back(node.at(k));
            dfs_lift(node.at(k));
        }
    }
}
