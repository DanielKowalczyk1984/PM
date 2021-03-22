#include "PricerSolverSimpleDP.hpp"
#include <fmt/core.h>
#include <gurobi_c++.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <vector>
#include "Instance.h"
#include "PricerSolverBase.hpp"
#include "scheduleset.h"

/**
 * Pricersolver for the TI index formulation
 */
PricerSolverSimpleDp::PricerSolverSimpleDp(const Instance& instance)
    : PricerSolverBase(instance),
      Hmax(instance.H_max),
      size_graph{},
      A(Hmax + 1),
      F(Hmax + 1),
      backward_F(Hmax + 1),
      TI_x(convex_constr_id * (Hmax + 1), GRBVar()),
      take(convex_constr_id * (Hmax + 1), 0),
      lp_x(convex_constr_id * (Hmax + 1), 0.0),
      solution_x(convex_constr_id * (Hmax + 1), 0.0) {
    init_table();
}

void PricerSolverSimpleDp::init_table() {
    backward_graph = std::vector<std::vector<Job*>>(Hmax + 1);
    forward_graph = std::vector<std::vector<Job*>>(Hmax + 1);

    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < convex_constr_id + 1; i++) {
            auto  j = i - 1;
            auto* job = jobs[j].get();

            if (t >= job->processing_time) {
                forward_graph[t].push_back(job);
                size_graph++;
            }

            if (t + job->processing_time <= Hmax) {
                backward_graph[t].push_back(job);
            }
        }
    }

    fmt::print("Number of arcs in TI formulation = {}\n", size_graph);
}

void PricerSolverSimpleDp::evaluate_nodes(double*                 pi,
                                          [[maybe_unused]] int    UB,
                                          [[maybe_unused]] double LB) {
    forward_evaluator(pi);
    backward_evaluator(pi);
}

void PricerSolverSimpleDp::evaluate_nodes([[maybe_unused]] double* pi) {
    forward_evaluator(pi);
    backward_evaluator(pi);
}

void PricerSolverSimpleDp::reduce_cost_fixing(double* pi, int UB, double LB) {
    evaluate_nodes(pi, UB, LB);
    std::span aux_pi{pi, reformulation_model.size()};
    int       counter = 0;
    int       x = 0;

    for (int t = 0; t < Hmax + 1; t++) {
        auto it = forward_graph[t].begin();
        while (it != forward_graph[t].end()) {
            double result = F[t - (*it)->processing_time] + value_Fj(t, *it) -
                            aux_pi[(*it)->job] + backward_F[t];
            if (LB + result + (convex_rhs - 1) * F[Hmax] > UB + RC_FIXING) {
                size_graph--;
                it = forward_graph[t].erase(it);
            } else {
                it++;
            }
        }

        x += forward_graph[t].size();

        auto iter = backward_graph[t].begin();
        while (iter != backward_graph[t].end()) {
            double result =
                F[t] + value_Fj(t + (*iter)->processing_time, *iter) -
                aux_pi[(*iter)->job] + backward_F[t + (*iter)->processing_time];
            if (LB + result + (convex_rhs - 1) * F[Hmax] > UB + RC_FIXING) {
                take[(*iter)->job * (Hmax + 1) + t] = false;
                // iter = backward_graph[t].erase(iter);
            } else {
                take[(*iter)->job * (Hmax + 1) + t] = true;
            }
            iter++;
        }

        counter += backward_graph[t].size();
    }

    fmt::print("new size of TI formulation = {} {} {}", size_graph, counter, x);
}

void PricerSolverSimpleDp::build_mip() {
    try {
        fmt::print("Building Mip model for the TI formulation\n");

        /** Constructing variables */
        for (int t = 0; t < Hmax + 1; t++) {
            for (auto& it : backward_graph[t]) {
                double cost = value_Fj(t + it->processing_time, it);
                double ub = take[(it->job) * (Hmax + 1) + t] ? 1.0 : 0.0;
                TI_x[it->job * (Hmax + 1) + t] =
                    model.addVar(0.0, ub, cost, GRB_BINARY);
            }
        }

        model.update();

        /** Assignment variables */
        std::vector<GRBLinExpr> assignment(convex_constr_id, GRBLinExpr());
        std::vector<char>       sense(convex_constr_id, GRB_EQUAL);
        std::vector<double>     rhs(convex_constr_id, 1.0);

        for (int t = 0; t <= Hmax; t++) {
            for (auto& it : backward_graph[t]) {
                assignment[it->job] += TI_x[it->job * (Hmax + 1) + t];
            }
        }

        std::unique_ptr<GRBConstr> assignment_constrs(
            model.addConstrs(assignment.data(), sense.data(), rhs.data(),
                             nullptr, convex_constr_id));

        model.update();

        std::vector<GRBLinExpr> interval_constr(Hmax + 1, GRBLinExpr());
        std::vector<char>       interval_sense(Hmax + 1);
        std::vector<double>     interval_rhs(Hmax + 1);

        for (auto t = 0UL; t <= Hmax; t++) {
            auto add_constraint = false;
            for (auto& it : backward_graph[t]) {
                for (int s = std::max(0UL, t - it->processing_time); s <= t;
                     s++) {
                    if (std::find(backward_graph[s].begin(),
                                  backward_graph[s].end(),
                                  it) != backward_graph[s].end()) {
                        interval_constr[t] += TI_x[(it->job) * (Hmax + 1) + s];
                        add_constraint = true;
                    }
                }
            }

            if (add_constraint) {
                interval_sense[t] = GRB_LESS_EQUAL;
                interval_rhs[t] = convex_rhs;

                model.addConstr(interval_constr[t], interval_sense[t],
                                interval_rhs[t]);
            }
        }

        model.update();

    } catch (GRBException& e) {
        std::cerr << e.getMessage() << '\n';
    }

    for (int t = 0; t <= Hmax; t++) {
        for (auto& it : backward_graph[t]) {
            TI_x[it->job * (Hmax + 1) + t].set(
                GRB_DoubleAttr_Start, solution_x[(it->job) * (Hmax + 1) + t]);
        }
    }

    model.write("ti_" + problem_name + "_" + std::to_string(convex_rhs) +
                "correct.lp");
    model.optimize();
}

void PricerSolverSimpleDp::forward_evaluator(double* _pi) {
    /** Initialisation */
    std::span aux_pi{_pi, reformulation_model.size()};
    F[0] = aux_pi[convex_constr_id];
    A[0] = nullptr;

    for (int t = 1; t < Hmax + 1; t++) {
        F[t] = DBL_MAX / 2;
        A[t] = nullptr;
    }

    /** Recursion */
    for (int t = 1; t < Hmax + 1; t++) {
        for (auto& it : forward_graph[t]) {
            if (F[t - it->processing_time] +
                    static_cast<double>(value_Fj(t, it)) - aux_pi[it->job] <=
                F[t]) {
                F[t] = F[t - it->processing_time] + value_Fj(t, it) -
                       aux_pi[it->job];
                A[t] = it;
            }
        }

        if (F[t - 1] <= F[t]) {
            F[t] = F[t - 1];
        }
    }
}

void PricerSolverSimpleDp::backward_evaluator(double* _pi) {
    backward_F[Hmax] = 0.0;
    std::span aux_pi{_pi, reformulation_model.size()};

    for (int t = 0; t < Hmax; t++) {
        backward_F[t] = DBL_MAX / 2;
    }

    for (auto t = Hmax - 1; t >= 0UL; t--) {
        for (auto& it : backward_graph[t]) {
            auto tt = t + it->processing_time;
            if (backward_F[tt] + static_cast<double>(value_Fj(tt, it)) -
                    aux_pi[it->job] <=
                backward_F[t]) {
                backward_F[t] = backward_F[tt] +
                                static_cast<double>(value_Fj(tt, it)) -
                                aux_pi[it->job];
            }
        }

        if (backward_F[t + 1] <= backward_F[t]) {
            backward_F[t] = backward_F[t + 1];
        }
    }
}

OptimalSolution<double> PricerSolverSimpleDp::pricing_algorithm(double* _pi) {
    OptimalSolution<double> opt_sol;
    opt_sol.cost = 0;
    std::vector<Job*> v;

    forward_evaluator(_pi);

    /** Find optimal solution */
    opt_sol.obj = std::numeric_limits<double>::max();

    for (int i = 0; i < Hmax + 1; i++) {
        if (F[i] < opt_sol.obj) {
            opt_sol.C_max = i;
            opt_sol.obj = F[i];
        }
    }

    auto t_min = opt_sol.C_max;

    /** Construct the solution */
    while (A[t_min] != nullptr) {
        Job* job = A[t_min];
        v.push_back(A[t_min]);
        opt_sol.cost += value_Fj(t_min, A[t_min]);
        t_min -= job->processing_time;
    }

    auto it = v.rbegin();

    for (; it != v.rend(); ++it) {
        opt_sol.jobs.push_back(*it);
    }

    /** Free the memory */
    return opt_sol;
}

OptimalSolution<double> PricerSolverSimpleDp::farkas_pricing(
    [[maybe_unused]] double* _pi) {
    OptimalSolution<double> opt_sol;

    return opt_sol;
}

void PricerSolverSimpleDp::construct_lp_sol_from_rmp(
    const double*                                    columns,
    const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets,
    int                                              num_columns) {
    std::span aux_cols{columns, schedule_sets.size()};
    // std::span aux_schedule_sets{schedule_sets->pdata, schedule_sets->len};
    std::fill(lp_x.begin(), lp_x.end(), 0.0);
    for (int k = 0; k < num_columns; k++) {
        if (aux_cols[k] > EPS_SOLVER) {
            auto* tmp = schedule_sets[k].get();
            int   t = 0;
            // std::span aux_jobs{tmp->job_list->pdata, tmp->job_list->len};
            for (auto& it : tmp->job_list) {
                lp_x[(it->job) * (Hmax + 1) + t] += aux_cols[k];
                t += it->processing_time;
            }
        }
    }
}

void PricerSolverSimpleDp::iterate_zdd() {}

void PricerSolverSimpleDp::create_dot_zdd([[maybe_unused]] const char* name) {
    boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>
        graph;
    for (int t = 0; t <= Hmax; t++) {
        boost::add_vertex(graph);
    }

    for (int t = 0; t < Hmax; t++) {
        for (auto& it : backward_graph[t]) {
            boost::add_edge(t, t + it->processing_time, graph);
        }
    }
    auto file_name = "TI_representation_" + problem_name + "_" +
                     std::to_string(convex_rhs) + ".gv";
    auto otf = std::ofstream(file_name);
    boost::write_graphviz(otf, graph);
    otf.close();
}

void PricerSolverSimpleDp::print_number_nodes_edges() {}

int PricerSolverSimpleDp::get_num_remove_nodes() {
    return 0;
}

int PricerSolverSimpleDp::get_num_remove_edges() {
    return 0;
}

size_t PricerSolverSimpleDp::get_nb_edges() {
    size_t nb_edges = 0u;
    for (int t = 0; t < Hmax + 1; t++) {
        nb_edges += forward_graph[t].size();
    }
    return nb_edges;
}

size_t PricerSolverSimpleDp::get_nb_vertices() {
    size_t nb_vertices = 0u;
    for (int t = 0; t < Hmax + 1; t++) {
        if (!forward_graph[t].empty()) {
            nb_vertices++;
        }
    }
    return nb_vertices;
}

int PricerSolverSimpleDp::get_num_layers() {
    return 0;
}

void PricerSolverSimpleDp::print_num_paths() {}

bool PricerSolverSimpleDp::check_schedule_set(
    [[maybe_unused]] const std::vector<Job*>& set) {
    // int t = 0;
    // for(unsigned int j = 0; j < set->len; j++) {
    //     Job* tmp_j = static_cast<Job*>() g_ptr_array_index()atic_cast<set,
    //     j>())); t += tmp_j->processing_time;

    //     if (t > Hmax) {
    //         return false;
    //     }

    //     auto it = std::find(forward_graph[t].begin(),
    //     forward_graph[t].end(),tmp_j); if (it == forward_graph[t].end()) {
    //         return false;
    //     }

    // }

    return true;
}
