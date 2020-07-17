#ifndef PRICER_SOLVER_ZDD_HPP
#define PRICER_SOLVER_ZDD_HPP

#include <NodeBddStructure.hpp>
#include "MipGraph.hpp"
#include "OptimalSolution.hpp"
#include "PricerSolverBase.hpp"
#include "wctprivate.h"

class PricerSolverZdd : public PricerSolverBase {
   public:
    std::unique_ptr<DdStructure<NodeZdd<double>>> decision_diagram;
    size_t                                        size_graph;

    int nb_removed_edges = 0;
    int nb_removed_nodes = 0;

    MipGraph                  mip_graph;
    std::unique_ptr<double[]> lp_x;
    std::unique_ptr<double[]> solution_x;

    PricerSolverZdd(GPtrArray* _jobs, int _num_machines,
                    GPtrArray* _ordered_jobs, const char* p_name);
    void         init_table() override;
    virtual void evaluate_nodes(double* pi, int UB, double LB) override = 0;

    void reduce_cost_fixing(double* pi, int UB, double LB) override;
    void remove_layers();
    void remove_edges();

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
    void   remove_layers_init();
    int*   get_take() override { return NULL; }

    OptimalSolution<double> farkas_pricing(double* pi) override;

    void add_constraint(Job* job, GPtrArray* list, int order) override;

    void update_constraints() override {}

    void update_reduced_costs_arcs(double* _pi, bool farkas = false) override {}

    void insert_constraints_lp(NodeData* pd) override {}

    void update_coeff_constraints() override {}

    // double compute_reduced_cost(const OptimalSolution<> &sol, double *pi,
    // double *lhs) override {
    //     double result = 0;

    //     return result;
    // }

    // double compute_lagrange(const OptimalSolution<> &sol, double *pi)
    // override {
    //     double result =0;

    //     return result;
    // }
};

#endif  // PRICER_SOLVER_ZDD_HPP
