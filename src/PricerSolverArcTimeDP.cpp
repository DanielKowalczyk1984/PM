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

    F = new double*[jobs->len + 1];
    for (unsigned i = 0; i < jobs->len + 1; ++i) {
        F[i] = new double[Hmax + 1]{};
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
    for (int t = 1; t < Hmax + 1; t++) {
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
        for (int t = 0; t <= Hmax - tmp->processing_time; ++t) {
            if (graph[j][t].empty()) {
                F[j][t] = DBL_MAX / 2;
            }
        }
    }
    std::cout << "Number of arcs in ATI formulation = " << size_graph << '\n';
}

void PricerSolverArcTimeDp::evaluate_nodes(double* pi, int UB, double LB) {
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
    }
    delete[] graph;

    for (int i = 0; i < n + 1; ++i) {
        delete[] F[i];
    }
    delete[] F;

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

OptimalSolution<double> PricerSolverArcTimeDp::pricing_algorithm(double* _pi) {
    OptimalSolution<double> sol(-_pi[n]);
    std::vector<Job*>       v;

    F[n][0] = _pi[n];
    double sigma = _pi[n];
    _pi[n] = 0;

    for (int t = 0; t < Hmax + 1; ++t) {
        for (int j = 0; j <= n; ++j) {
            Job* tmp = vector_jobs[j];
            A[j][t] = nullptr;
            B[j][t] = -1;
            F[j][t] = DBL_MAX / 2;
            job_iterator it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time) {
                F[j][t] = F[(*it)->job][t - p_matrix[(*it)->job][j]] +
                          value_Fj(t + tmp->processing_time, tmp) - _pi[j];
                A[j][t] = (*it);
                B[j][t] = t - p_matrix[(*it)->job][j];
                it++;
                while (it != graph[j][t].end()) {
                    double result = F[(*it)->job][t - p_matrix[(*it)->job][j]] +
                                    value_Fj(t + tmp->processing_time, tmp) -
                                    _pi[j];
                    if (F[j][t] >= result) {
                        F[j][t] = result;
                        A[j][t] = (*it);
                        B[j][t] = t - p_matrix[(*it)->job][j];
                    }
                    it++;
                }
            }
        }
    }

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
    _pi[n] = sigma;

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
    return 0u;
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