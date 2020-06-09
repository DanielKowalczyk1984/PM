#ifndef PRICER_SOLVER_BDD_HPP
#define PRICER_SOLVER_BDD_HPP
#include <NodeBddStructure.hpp>
#include "MipGraph.hpp"
#include "PricerSolverBase.hpp"
#include "ModelInterface.hpp"
#include <vector>

class PricerSolverBdd : public PricerSolverBase {
   public:
    std::unique_ptr<DdStructure<>> decision_diagram;
    ReformulationModel reformulation_model;
    size_t                         size_graph;

    int nb_removed_edges = 0;
    int nb_removed_nodes = 0;

    MipGraph                  mip_graph;
    std::unique_ptr<double[]> lp_x;
    std::unique_ptr<double[]> solution_x;
    int H_min;

    PricerSolverBdd(GPtrArray* _jobs, int _num_machines,
                    GPtrArray* _ordered_jobs, const char* p_name);
    PricerSolverBdd(GPtrArray* _jobs, int _num_machines,
                    GPtrArray* _ordered_jobs, int* _take_jobs, int _Hmax,
                    const char* _p_name);
    void         init_table() override;
    virtual void evaluate_nodes(double* pi, int UB, double LB) override = 0;
    void check_infeasible_arcs();
    void topdown_filtering();
    void bottum_up_filtering();
    void equivalent_paths_filtering();
    void print_representation_file();
    void calculate_H_min();
    void cleanup_arcs();

    void reduce_cost_fixing(double* pi, int UB, double LB) override;
    void remove_layers();
    void remove_edges();
    void remove_layers_init();

    void   construct_mipgraph();
    void   build_mip() override;
    void   construct_lp_sol_from_rmp(const double*    columns,
                                     const GPtrArray* schedule_sets,
                                     int              num_columns) override;
    void   represent_solution(Solution* sol) override;
    void   project_solution(Solution* sol) override;
    bool   check_schedule_set(GPtrArray* set) override;
    void   make_schedule_set_feasible(GPtrArray* set) override;
    void   disjunctive_inequality(double* x, Solution* sol) override;
    void   iterate_zdd() override;
    void   create_dot_zdd(const char* name) override;
    void   print_number_nodes_edges() override;
    int    get_num_remove_nodes() override;
    int    get_num_remove_edges() override;
    size_t get_nb_edges() override;
    size_t get_nb_vertices() override;
    int    get_num_layers() override;
    void   print_num_paths() override;
    void   add_constraint(Job* job, GPtrArray* list, int order) override;
    int* get_take() override {
        return NULL;
    };
    private:
    void add_inequality(std::vector<int> v1, std::vector<int> v2);
    void add_inequality(std::vector<int> v1);

    void update_constraints() override {

    }

    void update_reduced_costs_arcs() override {

    }
};
// int g_compare_duration(gconstpointer a, gconstpointer b);

#endif  // PRICER_SOLVER_BDD_HPP
