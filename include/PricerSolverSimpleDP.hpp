#ifndef PRICER_SOLVER_SIMPLE_DP_HPP
#define PRICER_SOLVER_SIMPLE_DP_HPP
#include <cstddef>               // for size_t
#include <memory>                // for make_unique, unique_ptr, shared_ptr
#include <vector>                // for vector
#include "Instance.h"            // for Instance
#include "PricerSolverBase.hpp"  // for PricerSolverBase
#include "PricingSolution.hpp"   // for PricingSolution
#include "gurobi_c++.h"          // for GRBVar
struct Job;
struct NodeData;
struct Column;
class PricerSolverSimpleDp : public PricerSolverBase {
   private:
    size_t                         Hmax;
    size_t                         size_graph;
    std::vector<Job*>              A;
    std::vector<double>            F;
    std::vector<double>            backward_F;
    std::vector<std::vector<Job*>> backward_graph;
    std::vector<std::vector<Job*>> forward_graph;
    std::vector<GRBVar>            TI_x;
    std::vector<int>               take;
    std::vector<double>            lp_x;
    std::vector<double>            solution_x;

   public:
    // PricerSolverSimpleDp(GPtrArray*  _jobs,
    //                      int         _num_machines,
    //                      int         _Hmax,
    //                      const char* p_name,
    //                      double      _UB);

    PricerSolverSimpleDp(const Instance& instance);
    PricerSolverSimpleDp& operator=(const PricerSolverSimpleDp&) = default;
    PricerSolverSimpleDp& operator=(PricerSolverSimpleDp&&) = default;
    PricerSolverSimpleDp(PricerSolverSimpleDp&&) = default;

    PricerSolverSimpleDp(const PricerSolverSimpleDp& src)
        : PricerSolverBase(src),
          Hmax(src.Hmax),
          size_graph(src.size_graph),
          A(Hmax + 1),
          F(Hmax + 1),
          backward_F(Hmax + 1),
          TI_x(convex_constr_id * (Hmax + 1), GRBVar()),
          take(convex_constr_id * (Hmax + 1)),
          lp_x(convex_constr_id * (Hmax + 1), 0.0),
          solution_x(convex_constr_id * (Hmax + 1)) {
        init_table();
    };

    ~PricerSolverSimpleDp() override = default;

    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return std::make_unique<PricerSolverSimpleDp>(*this);
    };

    void init_table();

    bool evaluate_nodes([[maybe_unused]] double* pi) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                               lambda,
        const std::vector<std::shared_ptr<Column>>& columns) override;

    size_t  get_nb_edges() override;
    size_t  get_nb_vertices() override;
    cpp_int print_num_paths() override;

    bool check_column(Column const* set) override;

    PricingSolution<double> pricing_algorithm(double* _pi) override;
    PricingSolution<double> farkas_pricing(double* _pi) override;
    void                    forward_evaluator(double* _pi);
    void                    backward_evaluator(double* _pi);

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
};

#endif  // PRICER_SOLVER_SIMPLE_DP_HPP
