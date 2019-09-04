#ifndef PRICER_SOLVER_BDD_HPP
#define PRICER_SOLVER_BDD_HPP
#include "PricerSolverBase.hpp"
#include "MipGraph.hpp"

class PricerSolverBdd : public PricerSolverBase
{
    public:
        std::unique_ptr<DdStructure<>> decision_diagram;
        size_t size_graph;

        int nb_removed_edges;
        int nb_removed_nodes;

        MipGraph mip_graph;
        std::unique_ptr<GRBEnv> env;
        std::unique_ptr<GRBModel> model;
        std::unique_ptr<double[]> lp_x;
        std::unique_ptr<double[]> solution_x;

        PricerSolverBdd(GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs);
        void init_table() override;
        virtual void evaluate_nodes(double* pi, int UB, double LB) override = 0;


        void reduce_cost_fixing(double* pi, int UB, double LB) override;
        void remove_layers();
        void remove_edges();
        void remove_layers_init();

        void construct_mipgraph();
        void build_mip() override;
        void construct_lp_sol_from_rmp(const double* columns, const GPtrArray* schedule_sets, int num_columns) override;
        void represent_solution(Solution* sol)  override;
        double* project_solution(Solution* sol) override;
        bool check_schedule_set(GPtrArray* set) override;
        void disjunctive_inequality(double* x, Solution* sol) override;
        void iterate_zdd() override;
        void create_dot_zdd(const char* name) override;
        void print_number_nodes_edges() override;
        int get_num_remove_nodes() override ;
        int get_num_remove_edges() override ;
        size_t get_size_data() override ;
        size_t get_size_graph() override ;
        int get_num_layers() override ;
        void print_num_paths() override;
        void add_constraint(Job* job, GPtrArray* list, int order) override;
};

#endif // PRICER_SOLVER_BDD_HPP
