#include "PricerSolverArcTimeDP.hpp"
#include <fmt/core.h>
#include <iostream>
#include <limits>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/take.hpp>
#include <vector>
#include "Instance.h"
#include "PricerSolverBase.hpp"
#include "gurobi_c++.h"
#include "gurobi_c.h"
#include "scheduleset.h"

PricerSolverArcTimeDp::PricerSolverArcTimeDp(const Instance& instance)
    : PricerSolverBase(instance),
      Hmax(instance.H_max),
      n(instance.nb_jobs),
      j0(),
      size_graph(0u),
      vector_jobs(),
      nb_edges_removed{},
      lp_x((n + 1) * (n + 1) * (Hmax + 1), 0.0),
      solution_x((n + 1) * (n + 1) * (Hmax + 1), 0.0) {
    for (auto i = 0UL; i < n; ++i) {
        vector_jobs.push_back((jobs)[i].get());
    }
    j0.job = n;
    vector_jobs.push_back(&j0);

    init_table();
}

#pragma clang diagnostic pop

void PricerSolverArcTimeDp::init_table() {
    graph = vector3d_jobs(n + 1);
    reversed_graph = vector3d_jobs(n + 1);
    for (auto j = 0UL; j < n + 1; j++) {
        reversed_graph[j] = vector2d_jobs(Hmax + 1);
    }

    forward_F = vector2d_dbl(jobs.size() + 1);
    backward_F = vector2d_dbl(jobs.size() + 1);
    for (unsigned i = 0; i < jobs.size() + 1; ++i) {
        forward_F[i] = vector1d_dbl(Hmax + 1, 0.0);
        backward_F[i] = vector1d_dbl(Hmax + 1, 0.0);
    }

    A = vector2d_jobs(jobs.size() + 1);
    for (unsigned i = 0; i < jobs.size() + 1; ++i) {
        A[i] = vector1d_jobs(Hmax + 1, nullptr);
    }

    B = vector2d_int(jobs.size() + 1);
    for (unsigned i = 0; i < jobs.size() + 1; ++i) {
        B[i] = vector1d_int(Hmax + 1);
    }

    arctime_x = vector_grb_var(jobs.size() + 1);
    for (auto i = 0UL; i < n + 1; i++) {
        arctime_x[i] = vector2d_grb_var(jobs.size() + 1);
        for (auto j = 0UL; j < n + 1; j++) {
            arctime_x[i][j] = vector1d_grb_var(Hmax + 1);
        }
    }

    for (auto j = 0UL; j < n; ++j) {
        graph[j] = vector2d_jobs(Hmax + 1);
        Job* tmp = jobs[j].get();
        for (int t = 0; t < Hmax + 1; t++) {
            for (auto& it : vector_jobs) {
                int p = (it->job != j) ? it->processing_time : 1;
                if (it != tmp && t - p >= 0 &&
                    t <= Hmax - tmp->processing_time) {
                    graph[j][t].push_back(it);
                    size_graph++;
                }
            }
        }
    }

    graph[n] = vector2d_jobs(Hmax + 1);
    for (int t = 0; t < Hmax + 1; t++) {
        for (auto& it : vector_jobs) {
            if (t >= it->processing_time) {
                graph[n][t].push_back(it);
                size_graph++;
            }
        }
    }

    /**
     * Remove all not needed arcs from the sets
     */
    for (auto i = 0UL; i < n - 1; ++i) {
        Job* tmp_i = vector_jobs[i];
        for (auto j = i + 1; j < n; ++j) {
            Job* tmp_j = vector_jobs[j];
            for (int t = tmp_i->processing_time;
                 t <= Hmax - tmp_j->processing_time; ++t) {
                if (delta1(i, j, t) >= 0) {
                    remove_arc(i, j, t);
                    size_graph--;
                } else {
                    remove_arc(
                        j, i,
                        t - tmp_i->processing_time + tmp_j->processing_time);
                    size_graph--;
                }
            }
        }
    }

    for (auto j = 0UL; j < n; ++j) {
        Job* tmp_j = vector_jobs[j];
        for (int t = tmp_j->processing_time; t < Hmax; ++t) {
            if (delta2(j, t) <= 0) {
                remove_arc(n, j, t - tmp_j->processing_time + 1);
                size_graph--;
            } else {
                remove_arc(j, n, t);
                size_graph--;
            }
        }
    }

    for (auto j = 0UL; j < n + 1; ++j) {
        Job* tmp = vector_jobs[j];
        for (int t = 0; t <= Hmax; ++t) {
            if (graph[j][t].empty()) {
                forward_F[j][t] = std::numeric_limits<double>::max() / 2;
            } else {
                for (auto& it : graph[j][t]) {
                    reversed_graph[it->job][t].push_back(tmp);
                }
            }
        }
    }

    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t <= Hmax; t++) {
            if (reversed_graph[j][t].empty()) {
                backward_F[j][t] = std::numeric_limits<double>::max() / 2;
            }
        }
    }

    fmt::print("Number of arcs in ATI formulation = {}\n", size_graph);
}

void PricerSolverArcTimeDp::evaluate_nodes(double*                 pi,
                                           [[maybe_unused]] int    UB,
                                           [[maybe_unused]] double LB) {
    forward_evaluator(pi);
    backward_evaluator(pi);
}

void PricerSolverArcTimeDp::evaluate_nodes([[maybe_unused]] double* pi) {
    forward_evaluator(pi);
    backward_evaluator(pi);
}

void PricerSolverArcTimeDp::build_mip() {
    fmt::print("Building Mip model for the arcTI formulation\n");

    /** Constructing variables */
    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
            for (auto& it : graph[j][t]) {
                double cost = vector_jobs[j]->weighted_tardiness_start(t);
                double UB = (it->job == vector_jobs[j]->job) ? convex_rhs : 1.0;
                auto   s =
                    (it->job == vector_jobs[j]->job) ? GRB_INTEGER : GRB_BINARY;
                arctime_x[it->job][j][t] = model.addVar(0.0, UB, cost, s);
            }
        }
    }

    model.update();

    /** Assignment variables */
    std::vector<GRBLinExpr> assignment(convex_constr_id, GRBLinExpr());
    std::vector<char>       sense(convex_constr_id, GRB_GREATER_EQUAL);
    std::vector<double>     rhs(convex_constr_id, 1.0);

    for (int j = 0; j < n; j++) {
        for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
            for (auto& it : graph[j][t]) {
                assignment[j] += arctime_x[it->job][j][t];
            }
        }
    }

    std::unique_ptr<GRBConstr> assignment_constrs(
        model.addConstrs(assignment.data(), sense.data(), rhs.data(), nullptr,
                         convex_constr_id));

    for (int i = 0; i < n; i++) {
        for (int t = 0; t <= Hmax - vector_jobs[i]->processing_time; t++) {
            GRBLinExpr expr{};
            for (auto& it : graph[i][t]) {
                expr += arctime_x[it->job][i][t];
            }

            for (auto& it :
                 reversed_graph[i][t + vector_jobs[i]->processing_time]) {
                expr -=
                    arctime_x[i][it->job][t + vector_jobs[i]->processing_time];
            }
            model.addConstr(expr, GRB_EQUAL, 0);
        }
    }

    for (int t = 0; t < Hmax; t++) {
        GRBLinExpr expr{};
        for (auto& it : graph[n][t]) {
            expr += arctime_x[it->job][n][t];
        }

        for (auto& it : reversed_graph[n][t + 1]) {
            expr -= arctime_x[n][it->job][t + 1];
        }

        model.addConstr(expr, GRB_EQUAL, 0);
    }

    GRBLinExpr expr{};
    for (auto& it : reversed_graph[n][0]) {
        expr += arctime_x[n][it->job][0];
    }
    model.addConstr(expr, GRB_EQUAL, convex_rhs);

    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
            for (auto& it : graph[j][t]) {
                arctime_x[it->job][j][t].set(
                    GRB_DoubleAttr_Start,
                    solution_x[(it->job) * (Hmax + 1) * (n + 1) +
                               j * (Hmax + 1) + t]);
                arctime_x[it->job][j][t].set(
                    GRB_DoubleAttr_PStart,
                    lp_x[(it->job) * (Hmax + 1) * (n + 1) + j * (Hmax + 1) +
                         t]);
            }
        }
    }

    model.write("ati_" + problem_name + "_" + std::to_string(convex_rhs) +
                ".lp");
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        for (int j = 0; j < n; j++) {
            for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
                for (auto& it : graph[j][t]) {
                    auto a = arctime_x[it->job][j][t].get(GRB_DoubleAttr_X);
                    if (a > 0) {
                        fmt::print("{} {} {}\n", j, t,
                                   jobs[j]->processing_time);
                    }
                }
            }
        }
    }
}

void PricerSolverArcTimeDp::reduce_cost_fixing(double* pi, int UB, double LB) {
    evaluate_nodes(pi, UB, LB);
    std::span aux_pi{pi, reformulation_model.size()};

    for (auto&& [j, tmp] :
         vector_jobs | ranges::views::take(n) | ranges::views::enumerate) {
        // Job* tmp = vector_jobs[j];
        for (int t = 0; t <= Hmax - tmp->processing_time; t++) {
            auto it = graph[j][t].begin();
            while (it != graph[j][t].end()) {
                double result =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    tmp->weighted_tardiness_start(t) - aux_pi[tmp->job] +
                    backward_F[tmp->job][t + tmp->processing_time];
                if (result + aux_pi[n] +
                        (convex_rhs - 1) * (forward_F[n][Hmax] + aux_pi[n]) +
                        LB >
                    UB - 1 + RC_FIXING) {
                    size_graph--;
                    nb_edges_removed++;
                    // Job* tmp_j = (*it);
                    auto pend = std::find(reversed_graph[(*it)->job][t].begin(),
                                          reversed_graph[(*it)->job][t].end(),
                                          vector_jobs[j]);
                    reversed_graph[(*it)->job][t].erase(pend);
                    graph[j][t].erase(it);
                } else {
                    it++;
                }
            }
        }
    }

    fmt::print(
        "size_graph after reduced cost fixing = {} and edges removed = {}\n",
        size_graph, nb_edges_removed);
}

PricerSolverArcTimeDp::~PricerSolverArcTimeDp() = default;

void PricerSolverArcTimeDp::backward_evaluator(double* _pi) {
    // backward_F[n][T] = 0;
    backward_F[n][Hmax] = 0.0;
    std::span aux_pi{_pi, reformulation_model.size()};

    for (int t = Hmax - 1; t >= 0; --t) {
        for (int i = 0; i <= n; ++i) {
            // Job* tmp = vector_jobs[i];
            backward_F[i][t] = ((i == n) && (t == Hmax))
                                   ? 0.0
                                   : std::numeric_limits<double>::max() / 2;
            auto it = reversed_graph[i][t].begin();

            if (!reversed_graph[i][t].empty() && t <= Hmax) {
                double reduced_cost = ((*it)->job == n)
                                          ? (*it)->weighted_tardiness_start(t)
                                          : (*it)->weighted_tardiness_start(t) -
                                                aux_pi[(*it)->job];
                int    tt = ((*it)->job != n) ? (*it)->processing_time
                            : (*it)->job == i ? 1
                                              : 0;
                backward_F[i][t] =
                    backward_F[(*it)->job][t + tt] + reduced_cost;
                it++;
                while (it != reversed_graph[i][t].end()) {
                    reduced_cost = ((*it)->job == n)
                                       ? (*it)->weighted_tardiness_start(t)
                                       : (*it)->weighted_tardiness_start(t) -
                                             aux_pi[(*it)->job];

                    tt = ((*it)->job != n) ? (*it)->processing_time
                         : (*it)->job == i ? 1
                                           : 0;
                    double result =
                        backward_F[(*it)->job][t + tt] + reduced_cost;

                    if (backward_F[i][t] >= result) {
                        backward_F[i][t] = result;
                    }
                    it++;
                }

            } else {
                backward_F[i][t] = std::numeric_limits<double>::max() / 2;
            }
        }
    }
}

void PricerSolverArcTimeDp::forward_evaluator(double* _pi) {
    forward_F[n][0] = 0;
    std::span aux_pi{_pi, reformulation_model.size()};

    for (int t = 0; t < Hmax + 1; ++t) {
        for (int j = 0; j <= n; ++j) {
            Job* tmp = vector_jobs[j];
            A[j][t] = nullptr;
            B[j][t] = -1;
            forward_F[j][t] = ((j == n) && (t == 0))
                                  ? 0
                                  : std::numeric_limits<double>::max() / 2;
            auto it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time) {
                double reduced_cost =
                    (j == n) ? tmp->weighted_tardiness_start(t)
                             : tmp->weighted_tardiness_start(t) - aux_pi[j];
                forward_F[j][t] =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    reduced_cost;

                A[j][t] = (*it);
                B[j][t] = t - (*it)->processing_time;
                it++;
                while (it != graph[j][t].end()) {
                    reduced_cost =
                        (j == n) ? tmp->weighted_tardiness_start(t)
                                 : tmp->weighted_tardiness_start(t) - aux_pi[j];
                    double result =
                        ((*it)->job != vector_jobs[j]->job)
                            ? forward_F[(*it)->job][t - (*it)->processing_time]
                            : forward_F[(*it)->job][t - 1];
                    result += reduced_cost;
                    if (forward_F[j][t] >= result) {
                        forward_F[j][t] = result;
                        A[j][t] = (*it);
                        B[j][t] = ((*it)->job != vector_jobs[j]->job)
                                      ? t - (*it)->processing_time
                                      : t - 1;
                    }
                    it++;
                }
            }
        }
    }
}

OptimalSolution<double> PricerSolverArcTimeDp::pricing_algorithm(double* _pi) {
    std::span               aux_pi{_pi, reformulation_model.size()};
    OptimalSolution<double> sol(aux_pi[n]);
    std::vector<Job*>       v;

    forward_evaluator(_pi);

    int job = n;
    int T = Hmax;

    while (T > 0) {
        int aux_job = A[job][T]->job;
        int aux_T = B[job][T];
        if (aux_job != n) {
            v.push_back(vector_jobs[aux_job]);
            sol.C_max += vector_jobs[aux_job]->processing_time;
            sol.cost += vector_jobs[aux_job]->weighted_tardiness_start(aux_T);
            sol.obj += vector_jobs[aux_job]->weighted_tardiness_start(aux_T) -
                       aux_pi[aux_job];
        }
        job = aux_job;
        T = aux_T;
    }

    sol.C_max = 0;
    for (auto it = v.rbegin(); it != v.rend(); it++) {
        sol.jobs.push_back(*it);
    }

    return sol;
}

OptimalSolution<double> PricerSolverArcTimeDp::farkas_pricing(
    [[maybe_unused]] double* _pi) {
    OptimalSolution<double> opt_sol;

    return opt_sol;
}

void PricerSolverArcTimeDp::construct_lp_sol_from_rmp(
    const double*                                    columns,
    const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets) {
    std::fill(lp_x.begin(), lp_x.end(), 0.0);
    // std::span aux_schedule_sets{schedule_sets->pdata, schedule_sets->len};
    std::span aux_cols{columns, schedule_sets.size()};
    for (auto k = 0UL; k < schedule_sets.size(); k++) {
        if (aux_cols[k] > 0.0) {
            auto  counter = 0UL;
            auto* tmp = schedule_sets[k].get();
            // std::span aux_jobs{tmp->job_list->pdata, tmp->job_list->len};
            int i = n;
            int t = 0;
            while (t < Hmax + 1) {
                Job* tmp_j = nullptr;
                int  j = n;

                if (counter < tmp->job_list.size()) {
                    tmp_j = tmp->job_list[counter];
                    j = tmp_j->job;
                }

                lp_x[i * (n + 1) * (Hmax + 1) + j * (Hmax + 1) + t] +=
                    aux_cols[k];

                if (tmp_j == nullptr) {
                    i = n;
                    t += 1;
                } else {
                    i = j;
                    t += tmp_j->processing_time;
                    counter++;
                }
            }
        }
    }
}

void PricerSolverArcTimeDp::iterate_zdd() {}

void PricerSolverArcTimeDp::create_dot_zdd([[maybe_unused]] const char* name) {}

void PricerSolverArcTimeDp::print_number_nodes_edges() {}

int PricerSolverArcTimeDp::get_num_remove_nodes() {
    return 0;
}

int PricerSolverArcTimeDp::get_num_remove_edges() {
    return nb_edges_removed;
}

size_t PricerSolverArcTimeDp::get_nb_edges() {
    size_t nb_edges = 0u;
    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t < Hmax + 1; t++) {
            nb_edges += graph[j][t].size();
        }
    }
    return nb_edges;
}

size_t PricerSolverArcTimeDp::get_nb_vertices() {
    size_t nb_vertices = 0u;
    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t < Hmax + 1; t++) {
            if (!graph[j][t].empty()) {
                nb_vertices++;
            }
        }
    }
    return nb_vertices;
}

int PricerSolverArcTimeDp::get_num_layers() {
    return 0;
}

void PricerSolverArcTimeDp::print_num_paths() {}

bool PricerSolverArcTimeDp::check_schedule_set(const std::vector<Job*>& set) {
    // std::span aux_set{set->pdata, set->len};
    size_t counter = set.size() - 1;
    int    i = n;
    int    t = 0;

    while (t < Hmax + 1) {
        Job* tmp_j = nullptr;
        int  j = n;

        if (counter < set.size()) {
            j = set[counter]->job;
        }

        if (std::find(graph[j][t].begin(), graph[j][t].end(), vector_jobs[i]) ==
            graph[j][t].end()) {
            return true;
        }

        if (tmp_j == nullptr) {
            i = n;
            t += 1;
        } else {
            i = j;
            t += tmp_j->processing_time;
            counter--;
        }
    }

    return false;
}
