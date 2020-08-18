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
#include "gurobi_c.h"

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
    model->set(GRB_IntParam_PoolSolutions, 1024);

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
    model->set(GRB_DoubleParam_TimeLimit, 100);

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
        root_node.lp_key = sigma.size();
        sigma.push_back(model->addVar(0.0, 1.0, 0.0, 'B',
                                      fmt::format("sigma_{}", sigma.size())));
        node_ids.push_back(root);
        dfs(root);

        GRBLinExpr          expr = 0;
        std::vector<double> coeffs(nb_jobs, 1.0);
        double              m =
            nb_machines % 2 == 0 ? 0.0 : static_cast<double>(nb_machines);
        q = model->addVar(0.0, GRB_INFINITY, 0.0, 'I');

        expr.addTerms(coeffs.data(), jobs_var.data(), nb_jobs);
        expr += m * sigma[0] - m * sigma[terminal_key] - 2.0 * q;
        model->addConstr(expr, '=', 1.0);
        model->update();

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "", 40,
                   e.getMessage());
    }
}

void ZeroHalfCuts::generate_cuts() {
    try {
        model->optimize();
        auto status = model->get(GRB_IntAttr_Status);

        if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE ||
            status == GRB_UNBOUNDED) {
            std::cout << "The model cannot be solved "
                      << "because it is infeasible or unbounded" << std::endl;
            return;
        }
        if (status != GRB_OPTIMAL) {
            std::cout << "Optimization was stopped with status " << status
                      << std::endl;
            return;
        }

        // Print number of solutions stored
        auto nb_solutions = model->get(GRB_IntAttr_SolCount);
        fmt::print("Number of solutions found: {}\n", nb_solutions);

        for (auto i = 0; i < nb_solutions; i++) {
            model->set(GRB_IntParam_SolutionNumber, i);
            int j = 0;
            std::for_each(jobs_var.cbegin(), jobs_var.cend(),
                          [&j](const auto& it) {
                              auto x = it.get(GRB_DoubleAttr_Xn);
                              if (x > 1e-4) {
                                  fmt::print("{} ", j);
                              }
                              j++;
                          });
            fmt::print("{}\n", q.get(GRB_DoubleAttr_Xn));
            int k = 0;
            std::for_each(sigma.cbegin(), sigma.cend(), [&](const auto& it) {
                auto x = it.get(GRB_DoubleAttr_Xn);
                if (x > 1e-4) {
                    auto& node = table->node(node_ids[k]);
                    if (node_ids[k] > 1) {
                        fmt::print("{} {}\n", node.get_nb_job(),
                                   node.get_weight());
                    } else {
                        fmt::print("TERMINAL NODE\n");
                    }
                }
                k++;
            });
            fmt::print("{0:-^{1}}\n", "", 40);
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
                child.lp_key = sigma.size();
                if (node.branch[i] == 1) {
                    terminal_key = child.lp_key;
                }

                sigma.push_back(model->addVar(
                    0.0, 1.0, 0.0, 'B', fmt::format("sigma_{}", sigma.size())));
                node_ids.push_back(node.branch[i]);
                auto& s_source = sigma[node.lp_key];
                auto& s_head = sigma.back();
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
                auto& s_source = sigma[node.lp_key];
                auto& s_head = sigma[child.lp_key];
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

void ZeroHalfCuts::dfs_evaluate(const NodeId& v) {
    auto& node = table->node(v);
    node.lp_visited = true;

    for (int i = 0; i < 2; i++) {
        if (node.lp_x[i] > 1e-6) {
            auto& child = table->node(node.branch[i]);
            if (!child.lp_visited) {
                child.lp_key = sigma.size();
                if (node.branch[i] == 1) {
                    terminal_key = child.lp_key;
                }

                sigma.push_back(model->addVar(
                    0.0, 1.0, 0.0, 'B', fmt::format("sigma_{}", sigma.size())));
                node_ids.push_back(node.branch[i]);
                auto& s_source = sigma[node.lp_key];
                auto& s_head = sigma.back();
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
                auto& s_source = sigma[node.lp_key];
                auto& s_head = sigma[child.lp_key];
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
