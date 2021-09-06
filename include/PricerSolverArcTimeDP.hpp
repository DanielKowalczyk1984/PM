#ifndef __PRICERSOLVERARCTIMEDP_H__
#define __PRICERSOLVERARCTIMEDP_H__

#include <algorithm>             // for remove
#include <cstddef>               // for size_t
#include <ext/alloc_traits.h>    // for __alloc_traits<>::value_type
#include <memory>                // for allocator_traits<>::value_type, make...
#include <vector>                // for vector
#include "Instance.h"            // for Instance
#include "Job.h"                 // for Job
#include "PricerSolverBase.hpp"  // for PricerSolverBase
#include "PricingSolution.hpp"   // for PricingSolution
#include "gurobi_c++.h"          // for GRBVar
struct NodeData;
struct Column;

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
    size_t         Hmax;
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

    bool evaluate_nodes([[maybe_unused]] double* pi) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                               lambda,
        const std::vector<std::shared_ptr<Column>>& columns) override;
    // void add_constraint(Job* job, GPtrArray* list, int order) override;

    PricingSolution<double> pricing_algorithm(double* _pi) override;
    PricingSolution<double> farkas_pricing(double* pi) override;

    size_t  get_nb_edges() override;
    size_t  get_nb_vertices() override;
    cpp_int print_num_paths() override;
    bool    check_column(Column const* set) override;

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

    void remove_arc(const size_t& i, const size_t& j, const size_t& t) {
        Job* tmp_i = vector_jobs[i];
        // auto it = graph[j][t].find(tmp_i);
        auto pend = std::remove(graph[j][t].begin(), graph[j][t].end(), tmp_i);
        graph[j][t].erase(pend);
    }

    int delta2(const size_t& j, const size_t& t) {
        Job* tmp_j = vector_jobs[j];
        return tmp_j->weighted_tardiness(t) - tmp_j->weighted_tardiness(t + 1);
    }

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
};
#endif // __PRICERSOLVERARCTIMEDP_H__