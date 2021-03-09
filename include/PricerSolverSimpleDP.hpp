#ifndef PRICER_SOLVER_SIMPLE_DP_HPP
#define PRICER_SOLVER_SIMPLE_DP_HPP
#include <memory>
#include "Instance.h"
#include "PricerSolverBase.hpp"
#include "gurobi_c++.h"

class PricerSolverSimpleDp : public PricerSolverBase {
   private:
    int                            Hmax;
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
    ~PricerSolverSimpleDp() override;

    std::unique_ptr<PricerSolverBase> clone() override { return nullptr; };

    void init_table();

    void evaluate_nodes([[maybe_unused]] double* pi,
                        [[maybe_unused]] int     UB,
                        [[maybe_unused]] double  LB) override;
    void evaluate_nodes([[maybe_unused]] double* pi) override;
    void reduce_cost_fixing([[maybe_unused]] double* pi,
                            [[maybe_unused]] int     UB,
                            [[maybe_unused]] double  LB) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                                    columns,
        const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets,
        int                                              num_columns) override;

    void   add_constraint(Job* job, GPtrArray* list, int order) override;
    void   iterate_zdd() override;
    void   create_dot_zdd(const char* name) override;
    void   print_number_nodes_edges() override;
    int    get_num_remove_nodes() override;
    int    get_num_remove_edges() override;
    size_t get_nb_edges() override;
    size_t get_nb_vertices() override;
    int    get_num_layers() override;
    void   print_num_paths() override;

    bool check_schedule_set(GPtrArray* set) override;
    bool check_schedule_set(const std::vector<Job*>& set) override;
    void make_schedule_set_feasible(GPtrArray* set) override;

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    forward_evaluator(double* _pi);
    void                    backward_evaluator(double* _pi);

    int* get_take() override {
        // int* tmp = take.data();
        // take = nullptr;
        return take.data();
    }

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
};

#endif  // PRICER_SOLVER_SIMPLE_DP_HPP
