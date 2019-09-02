#include "PricerSolverSimpleDP.hpp"

/**
 * Pricersolver for the TI index formulation
 */
PricerSolverSimpleDp::PricerSolverSimpleDp(GPtrArray* _jobs, int _num_machines,
                                           int _Hmax)
    : PricerSolverBase(_jobs, _num_machines),
      Hmax(_Hmax),
      size_graph(0u),
      A(new Job*[Hmax + 1]),
      F(new double[Hmax + 1]) {
    init_table();
}

void PricerSolverSimpleDp::init_table() {
    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < njobs + 1; i++) {
            int  j = i - 1;
            Job* job = reinterpret_cast<Job*>(g_ptr_array_index(jobs, j));

            if (t >= job->processing_time) {
                size_graph++;
            }
        }
    }

    std::cout << "Number of arcs in TI formulation = " << size_graph << '\n';
}

PricerSolverSimpleDp::~PricerSolverSimpleDp() {}

void PricerSolverSimpleDp::evaluate_nodes(double* pi, int UB, double LB) {
    return;
}

void PricerSolverSimpleDp::reduce_cost_fixing(double* pi, int UB, double LB) {
    return;
}

void PricerSolverSimpleDp::build_mip() {
    return;
}

OptimalSolution<double> PricerSolverSimpleDp::pricing_algorithm(double* _pi) {
    OptimalSolution<double> opt_sol;
    opt_sol.cost = 0;
    int               t_min = 0;
    std::vector<Job*> v;

    /** Initialisation */
    F[0] = -_pi[njobs];
    A[0] = nullptr;

    for (int t = 1; t < Hmax + 1; t++) {
        F[t] = -DBL_MAX / 2;
        A[t] = nullptr;
    }

    /** Recursion */
    for (int t = 0; t < Hmax + 1; t++) {
        for (int i = 1; i < njobs + 1; i++) {
            int  j = i - 1;
            Job* job = reinterpret_cast<Job*>(g_ptr_array_index(jobs, j));

            if (t >= job->processing_time) {
                if (F[t - job->processing_time] -
                        static_cast<double>(value_Fj(t, job)) + _pi[job->job] >=
                    F[t]) {
                    F[t] = F[t - job->processing_time] - value_Fj(t, job) +
                           _pi[job->job];
                    A[t] = job;
                }
            }
        }
    }

    /** Find optimal solution */
    opt_sol.obj = -DBL_MAX;

    for (int i = 0; i < Hmax + 1; i++) {
        if (F[i] > opt_sol.obj) {
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

void PricerSolverSimpleDp::construct_lp_sol_from_rmp(
    const double* columns, const GPtrArray* schedule_sets, int num_columns,
    double* x) {}

double* PricerSolverSimpleDp::project_solution(Solution* sol) {
    double* x = nullptr;

    return x;
}

void PricerSolverSimpleDp::represent_solution(Solution* sol) {}

void PricerSolverSimpleDp::add_constraint(Job* job, GPtrArray* list,
                                          int order) {}

void PricerSolverSimpleDp::iterate_zdd() {}

void PricerSolverSimpleDp::create_dot_zdd(const char* name) {}

void PricerSolverSimpleDp::print_number_nodes_edges() {}

int PricerSolverSimpleDp::get_num_remove_nodes() {
    return 0;
}

int PricerSolverSimpleDp::get_num_remove_edges() {
    return 0;
}

size_t PricerSolverSimpleDp::get_datasize() {
    return 0u;
}

size_t PricerSolverSimpleDp::get_size_graph() {
    return size_graph;
}

int PricerSolverSimpleDp::get_num_layers() {
    return 0;
}

void PricerSolverSimpleDp::print_num_paths() {
    // cout << "Number of paths: " <<
    // decision_diagram->evaluate(tdzdd::ZddCardinality<>()) << "\n";
}


bool PricerSolverSimpleDp::check_schedule_set(GPtrArray* set) {
    return true;
}

void PricerSolverSimpleDp::disjunctive_inequality(double* x, Solution* sol) {}