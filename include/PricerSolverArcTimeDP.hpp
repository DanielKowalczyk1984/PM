#include <cstddef>
#include <memory>
#include <vector>
#include "Instance.h"
#include "PricerSolverBase.hpp"
#include "gurobi_c++.h"

// using std::vector;
class PricerSolverArcTimeDp : public PricerSolverBase {
    using vector3d_jobs = std::vector<std::vector<std::vector<Job*>>>;
    using vector2d_jobs = std::vector<std::vector<Job*>>;
    using vector1d_jobs = std::vector<Job*>;
    using vector2d_dbl = std::vector<std::vector<double>>;
    using vector1d_dbl = std::vector<double>;
    using vector_grb_var = std::vector<std::vector<std::vector<GRBVar>>>;
    using vector2d_grb_var = std::vector<std::vector<GRBVar>>;
    using vector1d_grb_var = std::vector<GRBVar>;
    using vector2d_int = std::vector<std::vector<int>>;
    using vector1d_int = std::vector<int>;

   private:
    int            Hmax;
    size_t         n;
    size_t         size_graph;
    vector3d_jobs  graph;
    vector3d_jobs  reversed_graph;
    vector1d_jobs  vector_jobs;
    Job            j0;
    vector2d_dbl   forward_F;
    vector2d_dbl   backward_F;
    vector2d_jobs  A;
    vector2d_int   B;
    vector_grb_var arctime_x;
    int            nb_edges_removed;
    vector1d_dbl   lp_x;
    vector1d_dbl   solution_x;

   public:
    // PricerSolverArcTimeDp(GPtrArray*  _jobs,
    //                       int         _num_machines,
    //                       int         _Hmax,
    //                       const char* p_name,
    //                       double      _UB);

    PricerSolverArcTimeDp(const Instance& instance);
    PricerSolverArcTimeDp(const PricerSolverArcTimeDp&) = default;
    PricerSolverArcTimeDp(PricerSolverArcTimeDp&&) = default;
    PricerSolverArcTimeDp& operator=(const PricerSolverArcTimeDp&) = default;
    PricerSolverArcTimeDp& operator=(PricerSolverArcTimeDp&&) = default;

    ~PricerSolverArcTimeDp() override;
    // PricerSolverArcTimeDp(const PricerSolverArcTimeDp& src)
    //     : PricerSolverBase(src),
    //       Hmax(src.Hmax),
    //       n(src.n),
    //       j0(),
    //       vector_jobs(),
    //       nb_edges_removed(src.nb_edges_removed),
    //       lp_x((n + 1) * (n + 1) * (Hmax + 1), 0.0),
    //       solution_x((n + 1) * (n + 1) * (Hmax + 1), 0.0) {
    //     for (auto i = 0; i < n; ++i) {
    //         vector_jobs.push_back((*jobs)[i].get());
    //     }

    //     j0.job = n;
    //     vector_jobs.push_back(&j0);

    //     init_table();
    // }

    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverArcTimeDp>(*this);
    };

    void init_table();

    void evaluate_nodes([[maybe_unused]] double* pi) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                                    columns,
        const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets)
        override;
    // void add_constraint(Job* job, GPtrArray* list, int order) override;

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
    bool   check_schedule_set(const std::vector<Job*>& set) override;

    void forward_evaluator(double* pi);
    void backward_evaluator(double* _pi);

    int delta1(const size_t& i, const size_t& j, const int& t) {
        Job* tmp_i = vector_jobs[i];
        Job* tmp_j = vector_jobs[j];
        return (tmp_i->weighted_tardiness(t) +
                tmp_j->weighted_tardiness_start(t)) -
               (tmp_j->weighted_tardiness_start(t - tmp_i->processing_time) +
                tmp_i->weighted_tardiness(t + tmp_j->processing_time));
    }

    void remove_arc(const size_t& i, const size_t& j, const int& t) {
        Job* tmp_i = vector_jobs[i];
        // auto it = graph[j][t].find(tmp_i);
        auto pend = std::remove(graph[j][t].begin(), graph[j][t].end(), tmp_i);
        graph[j][t].erase(pend);
    }

    int delta2(const size_t& j, const int& t) {
        Job* tmp_j = vector_jobs[j];
        return tmp_j->weighted_tardiness(t) - tmp_j->weighted_tardiness(t + 1);
    }

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