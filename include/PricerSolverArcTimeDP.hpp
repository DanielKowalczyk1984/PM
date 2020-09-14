#include "PricerSolverBase.hpp"

class PricerSolverArcTimeDp : public PricerSolverBase {
   private:
    int                 Hmax;
    int                 n;
    size_t              size_graph;
    std::vector<Job*>** graph;
    std::vector<Job*>** reversed_graph;
    std::vector<Job*>   vector_jobs;
    Job                 j0;
    double**            forward_F;
    double**            backward_F;
    Job***              A;
    int**               B;
    GRBVar***           arctime_x;
    int                 nb_edges_removed;
    double*             lp_x;
    double*             solution_x;

   public:
    PricerSolverArcTimeDp(GPtrArray*  _jobs,
                          int         _num_machines,
                          int         _Hmax,
                          const char* p_name,
                          double      _UB);
    ~PricerSolverArcTimeDp();
    void init_table() override;

    void reduce_cost_fixing([[maybe_unused]] double* pi,
                            [[maybe_unused]] int     UB,
                            [[maybe_unused]] double  LB) override;
    void evaluate_nodes([[maybe_unused]] double* pi,
                        [[maybe_unused]] int     UB,
                        [[maybe_unused]] double  LB) override;

    void evaluate_nodes([[maybe_unused]] double* pi) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(const double*    columns,
                                   const GPtrArray* schedule_sets,
                                   int              num_columns) override;
    void represent_solution(Solution* sol) override;
    void project_solution(Solution* sol) override;
    void add_constraint(Job* job, GPtrArray* list, int order) override;

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* pi) override;

    void   iterate_zdd() override;
    void   create_dot_zdd(const char* name) override;
    void   print_number_nodes_edges() override;
    int    get_num_remove_nodes() override;
    int    get_num_remove_edges() override;
    size_t get_nb_edges() override;
    size_t get_nb_vertices() override;
    int    get_num_layers() override;
    void   print_num_paths() override;
    bool   check_schedule_set(GPtrArray* set) override;
    void   make_schedule_set_feasible(GPtrArray* set) override;

    void forward_evaluator(double* pi);
    void backward_evaluator(double* _pi);

    int delta1(const int& i, const int& j, const int& t) {
        Job* tmp_i = vector_jobs[i];
        Job* tmp_j = vector_jobs[j];
        return (value_Fj(t, tmp_i) +
                value_Fj(t + tmp_j->processing_time, tmp_j)) -
               (value_Fj(t + tmp_j->processing_time - tmp_i->processing_time,
                         tmp_j) +
                value_Fj(t + tmp_j->processing_time, tmp_i));
    }

    void remove_arc(const int& i, const int& j, const int& t) {
        Job* tmp_i = vector_jobs[i];
        // auto it = graph[j][t].find(tmp_i);
        auto pend = std::remove(graph[j][t].begin(), graph[j][t].end(), tmp_i);
        graph[j][t].erase(pend);
    }

    int delta2(const int& j, const int& t) {
        Job* tmp_j = vector_jobs[j];
        return value_Fj(t, tmp_j) - value_Fj(t + 1, tmp_j);
    }

    int* get_take() override { return NULL; }

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
    // double compute_reduced_cost(const OptimalSolution<>&s, double *pi, double
    // *lhs) override {
    //     double result = 0.0;

    //     return result;
    // }

    // double compute_lagrange(const OptimalSolution<> &sol, double *pi)
    // override {
    //     double result = 0.0;

    //     return result;
    // }
};