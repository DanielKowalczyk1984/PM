#ifndef PRICER_SOLVER_ZDD_HPP
#define PRICER_SOLVER_ZDD_HPP

#include "PricerSolverBase.hpp"
// #include "PricerEvaluate.hpp"
#include "MipGraph.hpp"
// #include "ZddNode.hpp"


class PricerSolverZdd : public PricerSolverBase
{
    public:
        std::unique_ptr<DdStructure<NodeZdd<double>>> decision_diagram;
        size_t size_graph;

        int nb_removed_edges;
        int nb_removed_nodes;

        MipGraph g;
        std::unique_ptr<GRBEnv> env;
        std::unique_ptr<GRBModel> model;

        PricerSolverZdd(GPtrArray* _jobs, int _num_machines, GPtrArray* _ordered_jobs);
        void init_table() override;
        virtual void evaluate_nodes(double* pi, int UB, double LB) override = 0;


        void reduce_cost_fixing(double* pi, int UB, double LB) override;
        void remove_layers();
        void remove_edges();

        void construct_mipgraph();
        void build_mip() override;
        void construct_lp_sol_from_rmp(const double* columns, const GPtrArray* schedule_sets, int num_columns, double* x) override;
        void represent_solution(Solution* sol)  override;
        double* project_solution(Solution* sol) override;
        bool check_schedule_set(GPtrArray* set) override;
        void disjunctive_inequality(double* x, Solution* sol) override;
        void iterate_zdd() override;
        void create_dot_zdd(const char* name) override;
        void print_number_nodes_edges() override;
        int get_num_remove_nodes() override ;
        int get_num_remove_edges() override ;
        size_t get_datasize() override ;
        size_t get_size_graph() override ;
        int get_num_layers() override ;
        void print_num_paths() override;
        double get_cost_edge(int idx) override;
        void remove_layers_init();


        void add_constraint(Job* job, GPtrArray* list, int order) override;
};


#endif // PRICER_SOLVER_ZDD_HPP
