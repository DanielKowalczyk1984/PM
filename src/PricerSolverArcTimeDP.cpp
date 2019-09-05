#include "PricerSolverArcTimeDP.hpp"

PricerSolverArcTimeDp::PricerSolverArcTimeDp(GPtrArray* _jobs,
                                             int _num_machines, int _Hmax)
    : PricerSolverBase(_jobs, _num_machines),
      Hmax(_Hmax),
      n(_jobs->len),
      size_graph(0u),
      vector_jobs(),
      env(new GRBEnv()),
      model(new GRBModel(*env)),
      num_edges_removed{},
      lp_x(new double[(n + 1) *(n + 1) * (Hmax + 1)]{}),
      solution_x(new double[(n + 1)* (n + 1)*(Hmax + 1)]{}) {
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

    arctime_x = new GRBVar**[jobs->len + 1];
    for (int i = 0; i < n + 1; i++) {
        arctime_x[i] = new GRBVar*[jobs->len + 1];
        for (int j = 0; j < n + 1; j++) {
            arctime_x[i][j] = new GRBVar[Hmax + 1];
        }
    }

    for (int j = 0; j < n; ++j) {
        graph[j] = new std::vector<Job*>[Hmax + 1];
        Job* tmp = reinterpret_cast<Job*>(g_ptr_array_index(jobs, j));
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

void PricerSolverArcTimeDp::evaluate_nodes(double* pi, int UB, double LB) {
    forward_evaluator(pi);
    backward_evaluator(pi);

    return;
}

void PricerSolverArcTimeDp::build_mip() {
    std::cout << "Building Mip model for the arcTI formulation\n";
    model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    // model->set(GRB_IntParam_Presolve, 2);
    // model->set(GRB_IntParam_VarBranch, 3);

    /** Constructing variables */
    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
            for (auto& it : graph[j][t]) {
                double cost = value_Fj(t + vector_jobs[j]->processing_time,
                                       vector_jobs[j]);
                double UB =
                    (it->job == vector_jobs[j]->job) ? num_machines : 1.0;
                auto s = (it->job == vector_jobs[j]->job) ? GRB_INTEGER : GRB_BINARY;
                arctime_x[it->job][j][t] =
                    model->addVar(0.0, UB, cost, s);
            }
        }
    }

    model->update();

    /** Assignment variables */
    std::unique_ptr<GRBLinExpr[]> assignment(new GRBLinExpr[nb_jobs]());
    std::unique_ptr<char[]>       sense(new char[nb_jobs]);
    std::unique_ptr<double[]>     rhs(new double[nb_jobs]);
    std::unique_ptr<double[]>     sum(new double[nb_jobs]);

    for (unsigned i = 0; i < jobs->len; ++i) {
        sense[i] = GRB_GREATER_EQUAL;
        rhs[i] = 1.0;
    }

    for (int j = 0; j < n; j++) {
        for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
            for (auto& it : graph[j][t]) {
                assignment[j] += arctime_x[it->job][j][t];
                if(t == 0 ) {
                    std::cout << (it->job) <<  " " << solution_x[(it->job)*(n + 1) * (Hmax + 1) + j* ( Hmax + 1) + t] << "\n";
                } 

                if( solution_x[(it->job)*(n + 1) * (Hmax + 1) + j* ( Hmax + 1) + t] > 0) {
                    std::cout << "OK OK OK " << (it->job) << " " << j << "\n";
                }
                sum[j] += solution_x[(it->job)*(n + 1) * (Hmax + 1) + j* ( Hmax + 1) + t];
            }
        }
    }

    for(int j = 0; j < n; j++) {
        std::cout << " sum = " << sum[j] << "\n";
    }

    std::unique_ptr<GRBConstr[]> assignment_constrs(model->addConstrs(
        assignment.get(), sense.get(), rhs.get(), nullptr, nb_jobs));

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
            model->addConstr(expr, GRB_EQUAL, 0);
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

        model->addConstr(expr, GRB_EQUAL, 0);
    }

    GRBLinExpr expr{};
    for (auto& it : reversed_graph[n][0]) {
        expr += arctime_x[n][it->job][0];
    }
    model->addConstr(expr, GRB_EQUAL, num_machines);
    
    for (int j = 0; j < n + 1; j++) {
        for (int t = 0; t <= Hmax - vector_jobs[j]->processing_time; t++) {
            for (auto& it : graph[j][t]) {
                arctime_x[it->job][j][t].set(GRB_DoubleAttr_Start, solution_x[(it->job)*(Hmax + 1)*(n + 1) + j*(Hmax + 1) + t]);
            }
        }
    }

    // model->optimize();
    return;
}

void PricerSolverArcTimeDp::reduce_cost_fixing(double* pi, int UB, double LB) {
    evaluate_nodes(pi, UB, LB);

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
                    num_edges_removed++;
                    // Job* tmp_j = (*it);
                    auto pend = std::find(reversed_graph[(*it)->job][t].begin(),reversed_graph[(*it)->job][t].end(),vector_jobs[j]);
                    reversed_graph[(*it)->job][t].erase(pend);
                    graph[j][t].erase(it);
                } else {
                    it++;
                }
            }
        }
    }

    std::cout << "size_graph after reduced cost fixing = " << size_graph
              << " and edges removed = " << num_edges_removed << "\n";

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

    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            delete[] arctime_x[i][j];
        }
        delete[] arctime_x[i];
    }
    delete[] arctime_x;
    delete[] lp_x;
    delete[] solution_x;
}

void PricerSolverArcTimeDp::backward_evaluator(double* _pi) {
    // backward_F[n][T] = 0;
    backward_F[n][Hmax] = 0;

    for (int t = Hmax - 1; t >= 0; --t) {
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
                int tt = ((*it)->job != n) ? (*it)->processing_time
                                           : (*it)->job == i ? 1 : 0;
                backward_F[i][t] =
                    backward_F[(*it)->job][t + tt] + reduced_cost;
                it++;
                while (it != reversed_graph[i][t].end()) {
                    reduced_cost =
                        ((*it)->job == n)
                            ? value_Fj(t + (*it)->processing_time, *it)
                            : value_Fj(t + (*it)->processing_time, *it) -
                                  _pi[(*it)->job];

                    tt = ((*it)->job != n) ? (*it)->processing_time
                                           : (*it)->job == i ? 1 : 0;
                    double result =
                        backward_F[(*it)->job][t + tt] + reduced_cost;

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
            auto it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time) {
                double reduced_cost =
                    (j == n) ? value_Fj(t + tmp->processing_time, tmp)
                             : value_Fj(t + tmp->processing_time, tmp) - _pi[j];
                forward_F[j][t] =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    reduced_cost;

                A[j][t] = (*it);
                B[j][t] = t - (*it)->processing_time;
                it++;
                while (it != graph[j][t].end()) {
                    reduced_cost =
                        (j == n)
                            ? value_Fj(t + tmp->processing_time, tmp)
                            : value_Fj(t + tmp->processing_time, tmp) - _pi[j];
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
    const double* columns, const GPtrArray* schedule_sets, int num_columns) {
    std::fill(lp_x, lp_x + get_size_data(), 0);
    for (int k = 0; k < num_columns; k++) {
        if (columns[k]) {
            size_t       counter = 0;
            ScheduleSet* tmp =
                (ScheduleSet*)g_ptr_array_index(schedule_sets, k);
            int i = n;
            int t = 0;
            while (t < Hmax + 1) {
                Job* tmp_j = nullptr;
                int  j = n;

                if (counter < tmp->job_list->len) {
                    tmp_j = (Job*)g_ptr_array_index(tmp->job_list, counter);
                    j = tmp_j->job;
                }

                lp_x[i * (n + 1) * (Hmax + 1) + j * (Hmax + 1) + t] += columns[k];

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

double* PricerSolverArcTimeDp::project_solution(Solution* sol) {
    // double* x = new double[(n + 1) * (n + 1) * (Hmax + 1)]{};
    std::fill(solution_x, solution_x + get_size_data(), 0.0);

    for (int it = 0; it < sol->nb_machines; it++) {
        GPtrArray* tmp = sol->part[it].machine;
        size_t     counter = tmp->len - 1;
        int        i = n;
        int        t = 0;
        while (t < Hmax + 1) {
            Job* tmp_j;
            int  j;
            if (counter < tmp->len) {
                tmp_j = (Job*)g_ptr_array_index(tmp, counter);
                j = tmp_j->job;
            } else {
                tmp_j = nullptr;
                j = n;
            }
            solution_x[i * (n + 1) * (Hmax + 1) + j * (Hmax + 1) + t] += 1.0;
std::cout << i << " test " << j << " " << t << " " << solution_x[i * (n + 1) * (Hmax + 1) + j * (Hmax + 1) + t] <<  "\n" ;

            if (tmp_j == nullptr) {
                i = n;
                t += 1;
            } else {
                i = j;
                t += tmp_j->processing_time;
                counter--;
            }
        }
        std::cout << "\n";
    }
getchar();
    return solution_x;
}

void PricerSolverArcTimeDp::represent_solution(Solution* sol) {
    double* x = project_solution(sol);

    for(int j = 0; j < n;j++) {
        for(int t = 0; t <= Hmax ; t++) {
            for (auto &it : graph[j][t]) {
                    std::cout <<solution_x[(it->job)*(n + 1)*(Hmax +1) + j*(Hmax + 1) + t] << " " << j << "\n";
            }
        }
    }
getchar();
    // delete[] x;
}

void PricerSolverArcTimeDp::add_constraint(Job* job, GPtrArray* list,
                                           int order) {}

void PricerSolverArcTimeDp::iterate_zdd() {}

void PricerSolverArcTimeDp::create_dot_zdd(const char* name) {}

void PricerSolverArcTimeDp::print_number_nodes_edges() {}

int PricerSolverArcTimeDp::get_num_remove_nodes() {
    return 0;
}

int PricerSolverArcTimeDp::get_num_remove_edges() {
    return num_edges_removed;
}

size_t PricerSolverArcTimeDp::get_size_data() {
    return (n + 1) * (n + 1) * (Hmax + 1);
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

bool PricerSolverArcTimeDp::check_schedule_set(GPtrArray* set) {
    size_t counter = set->len - 1;
    int    i = n;
    int    t = 0;

    while (t < Hmax + 1) {
        Job* tmp_j = nullptr;
        int  j = n;

        if (counter < set->len) {
            tmp_j = (Job*)g_ptr_array_index(set, counter);
            j = tmp_j->job;
        }

        if (std::find(graph[j][t].begin(), graph[j][t].end(), vector_jobs[i]) ==
            graph[j][t].end()) {
            return false;
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

    return true;
}

void PricerSolverArcTimeDp::disjunctive_inequality(double* x, Solution* sol) {}