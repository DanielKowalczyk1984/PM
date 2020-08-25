#ifndef PRICER_SOLVER_SIMPLE_DP_HPP
#define PRICER_SOLVER_SIMPLE_DP_HPP
#include <memory>
#include "PricerSolverBase.hpp"

class PricerSolverSimpleDp : public PricerSolverBase {
   private:
    int                       Hmax;
    size_t                    size_graph;
    std::unique_ptr<Job*[]>   A;
    std::unique_ptr<double[]> F;
    std::unique_ptr<double[]> backward_F;
    std::vector<Job*>*        backward_graph;
    std::vector<Job*>*        forward_graph;
    GRBVar*                   TI_x;
    int*                      take;
    double*                   lp_x;
    double*                   solution_x;

   public:
    PricerSolverSimpleDp(GPtrArray*  _jobs,
                         int         _num_machines,
                         int         _Hmax,
                         const char* p_name,
                         double      _UB);
    ~PricerSolverSimpleDp();
    void init_table() override;

    void evaluate_nodes([[maybe_unused]] double* pi,
                        [[maybe_unused]] int     UB,
                        [[maybe_unused]] double  LB) override;
    void reduce_cost_fixing([[maybe_unused]] double* pi,
                            [[maybe_unused]] int     UB,
                            [[maybe_unused]] double  LB) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(const double*    columns,
                                   const GPtrArray* schedule_sets,
                                   int              num_columns) override;
    void represent_solution(Solution* sol) override;
    void project_solution(Solution* sol) override;

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
    void make_schedule_set_feasible(GPtrArray* set) override;

    OptimalSolution<double> pricing_algorithm(double* _pi) override;
    OptimalSolution<double> farkas_pricing(double* _pi) override;
    void                    forward_evaluator(double* _pi);
    void                    backward_evaluator(double* _pi);

    int* get_take() override {
        int* tmp = take;
        take = NULL;
        return tmp;
    }

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
};

#endif  // PRICER_SOLVER_SIMPLE_DP_HPP
