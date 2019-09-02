#include "PricerSolverArcTimeDP.hpp"

PricerSolverArcTimeDp::PricerSolverArcTimeDp(GPtrArray* _jobs,
                                             int _num_machines, int _Hmax)
    : PricerSolverBase(_jobs, _num_machines),
      Hmax(_Hmax),
      n(_jobs->len),
      size_graph(0u),
      vector_jobs() {
    for (int i = 0; i < n; ++i) {
        vector_jobs.push_back(
            reinterpret_cast<Job*>(g_ptr_array_index(jobs, i)));
    }
    job_init(&j0, 0, 0, 0);
    j0.job = n;
    vector_jobs.push_back(&j0);

    init_table();
}

void PricerSolverArcTimeDp::init_table() {
    graph = new std::vector<Job*>*[n + 1];
    reversed_graph = new std::vector<Job*>*[n + 1];
    for (int j = 0; j < n + 1; j++) {
        reversed_graph[j] = new std::vector<Job*>[Hmax + 1];
    }

    forward_F = new double*[jobs->len + 1];
    backward_F = new double*[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        forward_F[i] = new double[Hmax + 1]{};
        backward_F[i] = new double[Hmax + 1]{};
    }

    A = new Job**[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        A[i] = new Job*[Hmax + 1];
    }

    B = new int*[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        B[i] = new int[Hmax + 1];
    }

    p_matrix = new int*[n + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        p_matrix[i] = new int[n + 1];
    }

    for (int i = 0; i < n; ++i) {
        int p = vector_jobs[i]->processing_time;
        for (int j = 0; j < n + 1; ++j) {
            p_matrix[i][j] = p;
        }
    }

    for (int j = 0; j < n + 1; ++j) {
        p_matrix[n][j] = (j == n) ? 1 : 0;
    }

    for (int j = 0; j < n; ++j) {
        graph[j] = new std::vector<Job*>[Hmax + 1];
        Job* tmp = reinterpret_cast<Job*>(g_ptr_array_index(jobs, j));
        for (int t = 0; t < Hmax + 1; t++) {
            for (auto& it : vector_jobs) {
                if (it != tmp && t - p_matrix[it->job][j] >= 0 &&
                    t <= Hmax - tmp->processing_time) {
                    graph[j][t].push_back(it);
                    size_graph++;
                }
            }
        }
    }

    graph[n] = new std::vector<Job*>[Hmax + 1];
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
    for (int i = 0; i < n - 1; ++i) {
        Job* tmp_i = vector_jobs[i];
        for (int j = i + 1; j < n; ++j) {
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

    for (int j = 0; j < n; ++j) {
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

    for (int j = 0; j < n + 1; ++j) {
        Job* tmp = vector_jobs[j];
        for (int t = 0; t <= Hmax; ++t) {
            if (graph[j][t].empty()) {
                forward_F[j][t] = DBL_MAX / 2;
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
                backward_F[j][t] = DBL_MAX / 2;
            }
        }
    }

    std::cout << "Number of arcs in ATI formulation = " << size_graph << '\n';
}

// bool PricerSolverArcTimeDp::constrainted_lagrangian(Job *job, int j, int t,
// double *pi, int UB, double LB) {
//     Job *tmp = vector_jobs[j];
//     double result = forward_F[job->job][t - job->processing_time] +
//     value_Fj(t + tmp->processing_time, tmp) + backward_F[tmp->job][t +
//     tmp->processing_time]; return result + pi[n] + (num_machines -
//     1)*(forward_F[n][Hmax] + pi[n]) + LB >= UB;
// }

void PricerSolverArcTimeDp::evaluate_nodes(double* pi, int UB, double LB) {
    forward_evaluator(pi);
    backward_evaluator(pi);

    for (int j = 0; j < n; j++) {
        Job* tmp = vector_jobs[j];
        for (int t = 0; t <= Hmax - tmp->processing_time; t++) {
            auto it = graph[j][t].begin();
            while (it != graph[j][t].end()) {
                double result =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    value_Fj(t + tmp->processing_time, tmp) - pi[tmp->job] +
                    backward_F[tmp->job][t + tmp->processing_time];
                if (result + pi[n] +
                        (num_machines - 1) * (forward_F[n][Hmax] + pi[n]) + LB >
                    UB - 1 + 0.00001) {
                    size_graph--;
                    graph[j][t].erase(it);
                } else {
                    it++;
                }
            }
        }
    }

    std::cout << "size_graph = " << size_graph << "\n";
    return;
}

void PricerSolverArcTimeDp::build_mip() {
    return;
}

void PricerSolverArcTimeDp::reduce_cost_fixing(double* pi, int UB, double LB) {
    return;
}

PricerSolverArcTimeDp::~PricerSolverArcTimeDp() {
    for (int i = 0; i < n + 1; ++i) {
        delete[] graph[i];
        delete[] reversed_graph[i];
    }
    delete[] graph;
    delete[] reversed_graph;
    for (int i = 0; i < n + 1; ++i) {
        delete[] forward_F[i];
        delete[] backward_F[i];
    }
    delete[] forward_F;
    delete[] backward_F;

    for (int i = 0; i < n + 1; ++i) {
        delete[] A[i];
    }
    delete[] A;

    for (int i = 0; i < n + 1; ++i) {
        delete[] B[i];
    }
    delete[] B;

    for (int i = 0; i < n + 1; ++i) {
        delete[] p_matrix[i];
    }
    delete[] p_matrix;
}

void PricerSolverArcTimeDp::backward_evaluator(double* _pi) {
    // backward_F[n][T] = 0;
    backward_F[n][Hmax] = 0;

    for (int t = Hmax; t >= 0; --t) {
        for (int i = 0; i <= n; ++i) {
            // Job* tmp = vector_jobs[i];
            backward_F[i][t] = ((i == n) && (t == Hmax)) ? 0.0 : DBL_MAX / 2;
            auto it = reversed_graph[i][t].begin();

            if (!reversed_graph[i][t].empty() && t <= Hmax) {
                double reduced_cost =
                    ((*it)->job == n)
                        ? value_Fj(t + (*it)->processing_time, *it)
                        : value_Fj(t + (*it)->processing_time, *it) -
                              _pi[(*it)->job];
                int tt = ((*it)->job != n) ? (*it)->processing_time : (*it)->job == i ? 1 : 0;
                backward_F[i][t] =
                    backward_F[(*it)->job][t + tt] +
                    reduced_cost;
                it++;
                while (it != reversed_graph[i][t].end()) {
                    reduced_cost =
                        ((*it)->job == n)
                            ? value_Fj(t + (*it)->processing_time, *it)
                            : value_Fj(t + (*it)->processing_time, *it) -
                                  _pi[(*it)->job];

                    tt = ((*it)->job != n) ? (*it)->processing_time : (*it)->job == i ? 1 : 0;
                    double result =
                        backward_F[(*it)->job][t + tt] +
                        reduced_cost;

                    if (backward_F[i][t] >= result) {
                        backward_F[i][t] = result;
                    }
                    it++;
                }

            } else {
                backward_F[i][t] = DBL_MAX / 2;
            }
        }
    }
}

void PricerSolverArcTimeDp::forward_evaluator(double* _pi) {
    forward_F[n][0] = 0;

    for (int t = 0; t < Hmax + 1; ++t) {
        for (int j = 0; j <= n; ++j) {
            Job* tmp = vector_jobs[j];
            A[j][t] = nullptr;
            B[j][t] = -1;
            forward_F[j][t] = ((j == n) && (t == 0)) ? 0 : DBL_MAX / 2;
            job_iterator it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time) {
                double reduced_cost =
                    (j == n) ? value_Fj(t + tmp->processing_time, tmp)
                             : value_Fj(t + tmp->processing_time, tmp) - _pi[j];
                forward_F[j][t] =
                    forward_F[(*it)->job][t - p_matrix[(*it)->job][j]] +
                    reduced_cost;

                A[j][t] = (*it);
                B[j][t] = t - p_matrix[(*it)->job][j];
                it++;
                while (it != graph[j][t].end()) {
                    reduced_cost =
                        (j == n)
                            ? value_Fj(t + tmp->processing_time, tmp)
                            : value_Fj(t + tmp->processing_time, tmp) - _pi[j];
                    double result =
                        forward_F[(*it)->job][t - p_matrix[(*it)->job][j]] +
                        reduced_cost;
                    if (forward_F[j][t] >= result) {
                        forward_F[j][t] = result;
                        A[j][t] = (*it);
                        B[j][t] = t - p_matrix[(*it)->job][j];
                    }
                    it++;
                }
            }
        }
    }
}

OptimalSolution<double> PricerSolverArcTimeDp::pricing_algorithm(double* _pi) {
    OptimalSolution<double> sol(-_pi[n]);
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
            sol.cost += value_Fj(aux_T + vector_jobs[aux_job]->processing_time,
                                 vector_jobs[aux_job]);
            sol.obj += _pi[aux_job] -
                       value_Fj(aux_T + vector_jobs[aux_job]->processing_time,
                                vector_jobs[aux_job]);
        }
        job = aux_job;
        T = aux_T;
    }

    sol.C_max = 0;
    for (auto& it : v) {
        g_ptr_array_add(sol.jobs, it);
    }

    return sol;
}

void PricerSolverArcTimeDp::construct_lp_sol_from_rmp(
    const double* columns, const GPtrArray* schedule_sets, int num_columns,
    double* x) {}

double* PricerSolverArcTimeDp::project_solution(Solution* sol) {
    double* x = nullptr;

    return x;
}

void PricerSolverArcTimeDp::represent_solution(Solution* sol) {}

void PricerSolverArcTimeDp::add_constraint(Job* job, GPtrArray* list,
                                           int order) {}

void PricerSolverArcTimeDp::iterate_zdd() {}

void PricerSolverArcTimeDp::create_dot_zdd(const char* name) {}

void PricerSolverArcTimeDp::print_number_nodes_edges() {}

int PricerSolverArcTimeDp::get_num_remove_nodes() {
    return 0;
}

int PricerSolverArcTimeDp::get_num_remove_edges() {
    return 0;
}

size_t PricerSolverArcTimeDp::get_datasize() {
    return 0u;
}

size_t PricerSolverArcTimeDp::get_size_graph() {
    return size_graph;
}

int PricerSolverArcTimeDp::get_num_layers() {
    return 0;
}

void PricerSolverArcTimeDp::print_num_paths() {
    // cout << "Number of paths: " <<
    // decision_diagram->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}

double PricerSolverArcTimeDp::get_cost_edge(int idx) {
    return 0.0;
}

bool PricerSolverArcTimeDp::check_schedule_set(GPtrArray* set) {
    return true;
}

void PricerSolverArcTimeDp::disjunctive_inequality(double* x, Solution* sol) {}