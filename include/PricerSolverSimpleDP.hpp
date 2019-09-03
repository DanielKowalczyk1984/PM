#ifndef PRICER_SOLVER_SIMPLE_DP_HPP
#define PRICER_SOLVER_SIMPLE_DP_HPP
#include "PricerSolverBase.hpp"

class PricerSolverSimpleDp : public PricerSolverBase {
private:
    int Hmax;
    size_t size_graph;
    std::unique_ptr<Job*[]> A;
    std::unique_ptr<double[]> F;
    
public:
    PricerSolverSimpleDp(GPtrArray *_jobs, int _num_machines, int _Hmax);
    ~PricerSolverSimpleDp();
    void init_table() override;
    
    void evaluate_nodes(double *pi, int UB, double LB) override;
    void reduce_cost_fixing(double *pi, int UB, double LB) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(const double *columns, const GPtrArray* schedule_sets, int num_columns, double *x) override;
    void represent_solution(Solution *sol)  override;
    double* project_solution(Solution *sol) override;

    void add_constraint(Job *job, GPtrArray *list, int order) override;
    void iterate_zdd() override;
    void create_dot_zdd(const char* name) override;
    void print_number_nodes_edges() override;
    int get_num_remove_nodes() override ;
    int get_num_remove_edges() override ;
    size_t get_size_data() override ;
    size_t get_size_graph() override ;
    int get_num_layers() override ;
    void print_num_paths() override;
    void disjunctive_inequality(double *x, Solution *sol) override;


    bool check_schedule_set(GPtrArray* set) override;

    OptimalSolution<double> pricing_algorithm(double *_pi) override;
};

#endif // PRICER_SOLVER_SIMPLE_DP_HPP


