#include "PricerSolverSimpleDP.hpp"
#include <scheduleset.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include "boost/graph/graphviz.hpp"

/**
 * Pricersolver for the TI index formulation
 */
PricerSolverSimpleDp::PricerSolverSimpleDp(GPtrArray*  _jobs,
                                           int         _num_machines,
                                           int         _Hmax,
                                           const char* p_name,
                                           double      _UB)
    : PricerSolverBase(_jobs, _num_machines, p_name, _UB),
      Hmax(_Hmax),
      size_graph(0u),
      A(new Job*[Hmax + 1]),
      F(new double[Hmax + 1]),
      backward_F(new double[Hmax + 1]),
      TI_x(new GRBVar[nb_jobs * (Hmax + 1)]),
      take(static_cast<int*>(malloc(nb_jobs * (Hmax + 1) * sizeof(int)))),
      lp_x(new double[nb_jobs * (Hmax + 1)]{}),
      solution_x(new double[nb_jobs * (Hmax + 1)]{}) {
    fill_int(take, nb_jobs * (Hmax + 1), 0);
    init_table();
}

void PricerSolverSimpleDp::init_table() {
    backward_graph = new std::vector<Job*>[Hmax + 1];
    forward_graph = new std::vector<Job*>[Hmax + 1];

    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < nb_jobs + 1; i++) {
            int  j = i - 1;
            Job* job = reinterpret_cast<Job*>(g_ptr_array_index(jobs, j));

            if (t >= job->processing_time) {
                forward_graph[t].push_back(job);
                size_graph++;
            }

            if (t + job->processing_time <= Hmax) {
                backward_graph[t].push_back(job);
            }
        }
    }

    std::cout << "Number of arcs in TI formulation = " << size_graph << '\n';
}

PricerSolverSimpleDp::~PricerSolverSimpleDp() {
    delete[] backward_graph;
    delete[] forward_graph;
    delete[] TI_x;
    if (take) {
        free(take);
    }
    delete[] lp_x;
    delete[] solution_x;
}

void PricerSolverSimpleDp::evaluate_nodes(double*                 pi,
                                          [[maybe_unused]] int    UB,
                                          [[maybe_unused]] double LB) {
    forward_evaluator(pi);
    backward_evaluator(pi);
    return;
}

void PricerSolverSimpleDp::reduce_cost_fixing(double* pi, int UB, double LB) {
    evaluate_nodes(pi, UB, LB);
    int counter = 0;
    int x = 0;

    for (int t = 0; t < Hmax + 1; t++) {
        auto it = forward_graph[t].begin();
        while (it != forward_graph[t].end()) {
            double result = F[t - (*it)->processing_time] + value_Fj(t, *it) -
                            pi[(*it)->job] + backward_F[t];
            if (LB + result + (num_machines - 1) * F[Hmax] > UB + 1e-6) {
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
                pi[(*iter)->job] + backward_F[t + (*iter)->processing_time];
            if (LB + result + (num_machines - 1) * F[Hmax] > UB + 1e-6) {
                take[(*iter)->job * (Hmax + 1) + t] = false;
                // iter = backward_graph[t].erase(iter);
            } else {
                take[(*iter)->job * (Hmax + 1) + t] = true;
            }
            iter++;
        }

        counter += backward_graph[t].size();
    }

    std::cout << "new size of TI formulation = " << size_graph << " " << counter
              << " " << x << "\n";
    return;
}

void PricerSolverSimpleDp::build_mip() {
    try {
        std::cout << "Building Mip model for the TI formulation\n";

        /** Constructing variables */
        for (int t = 0; t < Hmax + 1; t++) {
            for (auto& it : backward_graph[t]) {
                double cost = value_Fj(t + it->processing_time, it);
                double ub = take[(it->job) * (Hmax + 1) + t] ? 1.0 : 0.0;
                TI_x[it->job * (Hmax + 1) + t] =
                    model->addVar(0.0, ub, cost, GRB_BINARY);
            }
        }

        model->update();

        /** Assignment variables */
        std::unique_ptr<GRBLinExpr[]> assignment(new GRBLinExpr[nb_jobs]());
        std::unique_ptr<char[]>       sense(new char[nb_jobs]);
        std::unique_ptr<double[]>     rhs(new double[nb_jobs]);

        for (unsigned i = 0; i < jobs->len; ++i) {
            sense[i] = GRB_EQUAL;
            rhs[i] = 1.0;
        }

        for (int t = 0; t <= Hmax; t++) {
            for (auto& it : backward_graph[t]) {
                assignment[it->job] += TI_x[it->job * (Hmax + 1) + t];
            }
        }

        std::unique_ptr<GRBConstr[]> assignment_constrs(model->addConstrs(
            assignment.get(), sense.get(), rhs.get(), nullptr, nb_jobs));

        model->update();

        std::unique_ptr<GRBLinExpr[]> interval_constr(
            new GRBLinExpr[Hmax + 1]());
        std::unique_ptr<char[]>   interval_sense(new char[Hmax + 1]);
        std::unique_ptr<double[]> interval_rhs(new double[Hmax + 1]);

        for (int t = 0; t <= Hmax; t++) {
            auto add_constraint = false;
            for (auto& it : backward_graph[t]) {
                for (int s = std::max(0, t - it->processing_time); s <= t;
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
                interval_rhs[t] = num_machines;

                model->addConstr(interval_constr[t], interval_sense[t],
                                 interval_rhs[t]);
            }
        }

        model->update();

    } catch (GRBException e) {
        std::cerr << e.getMessage() << '\n';
    }

    for (int t = 0; t <= Hmax; t++) {
        for (auto& it : backward_graph[t]) {
            TI_x[it->job * (Hmax + 1) + t].set(
                GRB_DoubleAttr_Start, solution_x[(it->job) * (Hmax + 1) + t]);
        }
    }

    model->write("ti_" + problem_name + "_" + std::to_string(num_machines) +
                 "correct.lp");
    model->optimize();
    return;
}

void PricerSolverSimpleDp::forward_evaluator(double* _pi) {
    /** Initialisation */
    F[0] = _pi[nb_jobs];
    A[0] = nullptr;

    for (int t = 1; t < Hmax + 1; t++) {
        F[t] = DBL_MAX / 2;
        A[t] = nullptr;
    }

    /** Recursion */
    for (int t = 1; t < Hmax + 1; t++) {
        for (auto& it : forward_graph[t]) {
            if (F[t - it->processing_time] +
                    static_cast<double>(value_Fj(t, it)) - _pi[it->job] <=
                F[t]) {
                F[t] =
                    F[t - it->processing_time] + value_Fj(t, it) - _pi[it->job];
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

    for (int t = 0; t < Hmax; t++) {
        backward_F[t] = DBL_MAX / 2;
    }

    for (int t = Hmax - 1; t >= 0; t--) {
        for (auto& it : backward_graph[t]) {
            Job* job = it;
            int  tt = t + job->processing_time;
            if (backward_F[tt] + static_cast<double>(value_Fj(tt, job)) -
                    _pi[job->job] <=
                backward_F[t]) {
                backward_F[t] = backward_F[tt] +
                                static_cast<double>(value_Fj(tt, job)) -
                                _pi[job->job];
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
    int               t_min = 0;
    std::vector<Job*> v;

    forward_evaluator(_pi);

    /** Find optimal solution */
    opt_sol.obj = DBL_MAX;

    for (int i = 0; i < Hmax + 1; i++) {
        if (F[i] < opt_sol.obj) {
            opt_sol.C_max = i;
            opt_sol.obj = F[i];
        }
    }

    t_min = opt_sol.C_max;

    /** Construct the solution */
    while (A[t_min] != nullptr) {
        Job* job = A[t_min];
        v.push_back(A[t_min]);
        opt_sol.cost += value_Fj(t_min, A[t_min]);
        t_min -= job->processing_time;
    }

    std::vector<Job*>::reverse_iterator it = v.rbegin();

    for (; it != v.rend(); ++it) {
        g_ptr_array_add(opt_sol.jobs, *it);
    }

    /** Free the memory */
    return opt_sol;
}

OptimalSolution<double> PricerSolverSimpleDp::farkas_pricing([
    [maybe_unused]] double* _pi) {
    OptimalSolution<double> opt_sol;

    return opt_sol;
}

void PricerSolverSimpleDp::construct_lp_sol_from_rmp(
    const double*    columns,
    const GPtrArray* schedule_sets,
    int              num_columns) {
    std::fill(lp_x, lp_x + nb_jobs * (Hmax + 1), 0.0);
    for (int k = 0; k < num_columns; k++) {
        if (columns[k] > 0.00001) {
            ScheduleSet* tmp =
                (ScheduleSet*)g_ptr_array_index(schedule_sets, k);
            int t = 0;
            for (size_t l = 0; l < tmp->job_list->len; l++) {
                Job* tmp_j = (Job*)g_ptr_array_index(tmp->job_list, l);
                lp_x[(tmp_j->job) * (Hmax + 1) + t] += columns[k];
                t += tmp_j->processing_time;
            }
        }
    }
}

void PricerSolverSimpleDp::project_solution(Solution* sol) {
    std::fill(solution_x, solution_x + nb_jobs * (Hmax + 1), 0.0);

    for (int it = 0; it < sol->nb_machines; it++) {
        GPtrArray* tmp = sol->part[it].machine;
        int        t = 0;
        for (size_t l = 0; l < tmp->len; l++) {
            Job* tmp_j = (Job*)g_ptr_array_index(tmp, l);
            solution_x[(tmp_j->job) * (Hmax + 1) + t] += 1.0;
            t += tmp_j->processing_time;
        }
    }
}

void PricerSolverSimpleDp::represent_solution(Solution* sol) {
    project_solution(sol);
}

void PricerSolverSimpleDp::add_constraint([[maybe_unused]] Job*       job,
                                          [[maybe_unused]] GPtrArray* list,
                                          [[maybe_unused]] int        order) {}

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
                     std::to_string(num_machines) + ".gv";
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

bool PricerSolverSimpleDp::check_schedule_set([[maybe_unused]] GPtrArray* set) {
    // int t = 0;
    // for(unsigned int j = 0; j < set->len; j++) {
    //     Job* tmp_j = (Job*) g_ptr_array_index(set, j);
    //     t += tmp_j->processing_time;

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

void PricerSolverSimpleDp::make_schedule_set_feasible([
    [maybe_unused]] GPtrArray* set) {}
