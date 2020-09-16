#include "ZeroHalfCuts.hpp"
#include <algorithm>
#include <memory>
#include <vector>
#include "ModelInterface.hpp"
#include "NodeBddTable.hpp"
#include "NodeId.hpp"
#include "fmt/core.h"
#include "fmt/format.h"
#include "gurobi_c++.h"

ZeroHalfCuts::ZeroHalfCuts(int                 _nb_jobs,
                           int                 _nb_machines,
                           ReformulationModel* _rmp_model,
                           NodeId&             _root,
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
    model->set(GRB_IntParam_PoolSolutions, 124);

    // Limit the search space by setting a gap for the worst possible solution
    // that will be accepted
    model->set(GRB_DoubleParam_PoolGap, 0.10);

    // do a systematic search for the k-best solutions
    model->set(GRB_IntParam_PoolSearchMode, 2);
    model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->set(GRB_IntParam_Presolve, 2);
    model->set(GRB_DoubleParam_Heuristics, 0.5);
    model->set(GRB_IntParam_MIPFocus, 1);
    model->set(GRB_DoubleParam_TimeLimit, 10);

    generate_model();
}

bool ZeroHalfCuts::add_cuts() {
    return true;
}

void ZeroHalfCuts::generate_model() {
    try {
        for (int i = 0; i < nb_jobs; i++) {
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
        expr += m * root_node.sigma - m * terminal_node.sigma - 2.0 * q;
        model->addConstr(expr, '=', 1.0);
        model->update();
        init_table();

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "", 40,
                   e.getMessage());
    }
}

void ZeroHalfCuts::init_table() {
    for (auto i = root.row(); i >= 0; i--) {
        for (auto& it : (*table)[i]) {
            it.in_degree_0 = 0;
            it.in_degree_1 = 0;
            it.in_edges[0].clear();
            it.in_edges[1].clear();
            it.coeff_cut[0] = 0.0;
            it.coeff_cut[1] = 0.0;
        }
    }

    for (auto i = root.row(); i > 0; i--) {
        for (auto& it : (*table)[i]) {
            auto& n0 = table->node(it.branch[0]);
            auto& n1 = table->node(it.branch[1]);

            n0.in_edges[0].push_back(it.ptr_node_id);
            n0.in_degree_0++;
            n1.in_edges[1].push_back(it.ptr_node_id);
            n1.in_degree_1++;
        }
    }
}

void ZeroHalfCuts::init_coeff_cut() {
    for (auto& it : node_ids) {
        auto& node = table->node(it);
        init_coeff_node(node);
    }

    for (auto& iter : node_ids_lift) {
        auto& node = table->node(iter);
        init_coeff_node(node);
    }
    node_ids_lift.clear();
}

void ZeroHalfCuts::init_coeff_node(NodeBdd<>& node) {
    for (int k = 0; k < 2; k++) {
        node.coeff_cut[k] = 0.0;
        for (auto& it : node.in_edges[k]) {
            auto aux = it.lock();
            if (aux) {
                auto& aux_node = table->node(*aux);
                aux_node.coeff_cut[k] = 0.0;
            }
        }
    }
}

void ZeroHalfCuts::construct_cut() {
    GenericData* data = new GenericData();

    auto add_coeff_constr = [&](const auto& it) {
        auto& node = table->node(it);

        for (int k = 0; k < 2; k++) {
            auto coeff = floor(node.coeff_cut[k] / 2.0);
            if (coeff > 1e-4)
                data->add_coeff_hash_table(node.get_nb_job(), node.get_weight(),
                                           k, -coeff);
            for (auto& iter : node.in_edges[k]) {
                auto aux = iter.lock();
                if (aux) {
                    auto& aux_node = table->node(*aux);
                    auto  coeff_in = floor(aux_node.coeff_cut[k] / 2);
                    if (coeff_in < -1e-4) {
                        data->add_coeff_hash_table(aux_node.get_nb_job(),
                                                   aux_node.get_weight(), k,
                                                   -coeff_in);
                    }
                }
            }
        }
    };

    auto print_node_ids = [&](const auto& it) {
        auto& node = table->node(it);
        if (it > 1 && node.sigma.get(GRB_DoubleAttr_Xn) > 0.0) {
            fmt::print("node {} {} {} | ", node.get_nb_job(), node.get_weight(),
                       node.sigma.get(GRB_DoubleAttr_Xn));
        } else if (node.sigma.get(GRB_DoubleAttr_Xn) > 0.0) {
            fmt::print("Terminal Node | ");
        }
    };

    std::for_each(node_ids.begin(), node_ids.end(), add_coeff_constr);
    std::for_each(node_ids_lift.begin(), node_ids_lift.end(), add_coeff_constr);
    // std::for_each(node_ids.begin(), node_ids.end(), print_node_ids);
    fmt::print("\n");

    auto& root_node = table->node(root);
    auto& terminal_node = table->node(1);
    auto  rhs = (root_node.sigma.get(GRB_DoubleAttr_Xn) > 1e-6)
                   ? static_cast<double>(nb_machines)
                   : 0.0;
    rhs += (terminal_node.sigma.get(GRB_DoubleAttr_Xn) > 1e-6)
               ? -static_cast<double>(nb_machines)
               : 0.0;
    std::shared_ptr<ConstraintGeneric> constr{
        std::make_shared<ConstraintGeneric>(
            data, -floor((2.0 * q.get(GRB_DoubleAttr_Xn) + 1 + rhs) / 2.0))};
    data->list_coeff();
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
            fmt::print("Optimization was stopped with status {}", status);
        }

        // Print number of solutions stored
        auto nb_solutions = model->get(GRB_IntAttr_SolCount);
        fmt::print("Number of solutions found: {}\n", nb_solutions);

        for (auto i = 0; i < std::min(5, nb_solutions); i++) {
            model->set(GRB_IntParam_SolutionNumber, i);
            init_coeff_cut();

            auto calc_coeff_cut = [&](const auto& it) {
                auto& node = table->node(it);
                auto  x = node.sigma.get(GRB_DoubleAttr_Xn);
                if (x > 1e-4) {
                    if (it > 1) {
                        node.coeff_cut[0] += 1.0;
                        node.coeff_cut[1] += 1.0;
                        auto y =
                            jobs_var[node.get_nb_job()].get(GRB_DoubleAttr_Xn);
                        if (y > 1e-4) {
                            node.coeff_cut[1] += 1.0;
                        }
                    }

                    for (int j = 0; j < 2; j++) {
                        for (auto& it_aux : node.in_edges[j]) {
                            auto aux = it_aux.lock();
                            if (aux) {
                                auto& aux_node = table->node(*aux);
                                aux_node.coeff_cut[j] -= 1.0;
                                if (j) {
                                    auto y =
                                        jobs_var[aux_node.get_nb_job()].get(
                                            GRB_DoubleAttr_Xn);
                                    if (y == 1.0) {
                                        aux_node.coeff_cut[j] += 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            };

            std::for_each(node_ids.begin(), node_ids.end(), calc_coeff_cut);

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
            //             if (j) {
            //                 auto y = jobs_var[child.get_nb_job()].get(
            //                     GRB_DoubleAttr_Xn);
            //                 if (y > 1e-4) {
            //                     child.coeff_cut[j] += 1.0;
            //                 }
            //             }
            //             for (auto& it : child.in_edges[j]) {
            //                 auto aux = it.lock();
            //                 if (aux) {
            //                     auto& aux_node = table->node(*aux);
            //                     aux_node.coeff_cut[j] -= 1.0;
            //                     auto y = jobs_var[aux_node.get_nb_job()].get(
            //                         GRB_DoubleAttr_Xn);
            //                     if (y == 1.0 && !aux_node.lp_visited) {
            //                         aux_node.coeff_cut[j] += 1.0;
            //                     }
            //                 }
            //             }
            //         }
            //         node_ids_lift.push_back(node.branch[0]);
            //         dfs_lift(node.branch[0]);
            //     }
            // }
            construct_cut();
        }

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "", 40,
                   e.getMessage());
    }
}

void ZeroHalfCuts::dfs(const NodeId& v) {
    auto& node = table->node(v);
    node.lp_visited = true;

    for (int i = 0; i < 2; i++) {
        if (node.lp_x[i] > 1e-6) {
            auto& child = table->node(node.branch[i]);
            if (!child.lp_visited) {
                if (node.branch[i] == 1) {
                    child.sigma = model->addVar(0.0, 1.0, 0.0, 'B',
                                                fmt::format("sigma_terminal"));
                } else {
                    child.sigma = model->addVar(
                        0.0, 1.0, 0.0, 'B',
                        fmt::format("sigma_{}_{}", child.get_nb_job(),
                                    child.get_weight()));
                }

                node_ids.push_back(node.branch[i]);
                auto& s_source = node.sigma;
                auto& s_head = child.sigma;
                auto  str_y = fmt::format("y_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto  str_r = fmt::format("r_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto& y = node.y[i] =
                    model->addVar(0.0, GRB_INFINITY, node.lp_x[i], 'B', str_y);
                auto& r = node.r[i] =
                    model->addVar(0.0, GRB_INFINITY, 0.0, 'I', str_r);
                GRBLinExpr expr = s_head - s_source - y - 2.0 * r;
                if (i) {
                    expr += jobs_var[node.get_nb_job()];
                }
                model->addConstr(expr, '=', 0.0);
                dfs(node.branch[i]);
            } else {
                auto& s_source = node.sigma;
                auto& s_head = child.sigma;
                auto  str_y = fmt::format("y_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto  str_r = fmt::format("r_{}_{}", node.get_nb_job(),
                                         node.get_weight());
                auto& y = node.y[i] =
                    model->addVar(0.0, GRB_INFINITY, node.lp_x[i], 'B', str_y);
                auto& r = node.r[i] =
                    model->addVar(0.0, GRB_INFINITY, 0.0, 'I', str_r);
                GRBLinExpr expr = s_head - s_source - y - 2.0 * r;
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

    for (int k = 0; k < 2; k++) {
        if (node.branch[k] <= 0) {
            continue;
        }
        auto& child_node = table->node(node.branch[k]);

        if (node.lp_x[k] == 0.0 && !child_node.lp_visited) {
            child_node.coeff_cut[0] += 1.0;
            child_node.coeff_cut[1] += 1.0;
            auto y = jobs_var[child_node.get_nb_job()].get(GRB_DoubleAttr_Xn);
            if (y == 1.0) {
                child_node.coeff_cut[1] += 1.0;
            }

            for (int j = 0; j < 2; j++) {
                for (auto& it : child_node.in_edges[j]) {
                    auto aux = it.lock();
                    if (aux) {
                        auto& aux_node = table->node(*aux);
                        aux_node.coeff_cut[j] -= 1.0;
                        if (j) {
                            auto y_parrent =
                                jobs_var[aux_node.get_nb_job()].get(
                                    GRB_DoubleAttr_Xn);
                            if (y == 1.0) {
                                aux_node.coeff_cut[j] += 1.0;
                            }
                        }
                    }
                }
            }
            node_ids_lift.push_back(node.branch[k]);
            dfs_lift(node.branch[k]);
        }
    }
}
