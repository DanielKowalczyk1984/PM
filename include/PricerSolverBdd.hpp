#ifndef PRICER_SOLVER_BDD_HPP
#define PRICER_SOLVER_BDD_HPP
#include <NodeBddStructure.hpp>
#include <memory>
#include <unordered_map>
#include <vector>
#include "MipGraph.hpp"
#include "ModelInterface.hpp"
#include "OptimalSolution.hpp"
#include "PricerSolverBase.hpp"
#include "wctprivate.h"

class PricerSolverBdd : public PricerSolverBase {
    std::unique_ptr<DdStructure<>> decision_diagram;
    size_t                         size_graph;
    int                            nb_removed_edges{};
    int                            nb_removed_nodes{};

    GPtrArray* ordered_jobs;

    MipGraph mip_graph;

    std::vector<std::vector<std::weak_ptr<NodeId>>> node_ids;

    OriginalModel<> original_model;

    std::unordered_map<int, std::vector<std::weak_ptr<NodeId>>> t_in;
    std::unordered_map<int, std::vector<std::weak_ptr<NodeId>>> t_out;

    int H_min;
    int H_max;

   public:
    PricerSolverBdd(GPtrArray*  _jobs,
                    int         _num_machines,
                    GPtrArray*  _ordered_jobs,
                    const char* p_name,
                    int         _Hmax,
                    int*        _take_jobs,
                    double      _UB);

    PricerSolverBdd(const PricerSolverBdd& src)
        : PricerSolverBase(src),
          decision_diagram(new DdStructure<>(*src.decision_diagram)),
          size_graph(src.size_graph),
          nb_removed_edges(src.nb_removed_edges),
          nb_removed_nodes(src.nb_removed_nodes),
          ordered_jobs(src.ordered_jobs),
          mip_graph(src.mip_graph),
          original_model(src.original_model) {}

    virtual void evaluate_nodes(double* pi, int UB, double LB) override = 0;

    void check_infeasible_arcs();
    void topdown_filtering();
    void bottum_up_filtering();
    void equivalent_paths_filtering();
    void print_representation_file();
    void calculate_H_min();
    void cleanup_arcs();

    void remove_layers();
    void remove_edges();
    void remove_layers_init();
    void construct_mipgraph();
    void init_coeff_constraints();

    bool   check_schedule_set(GPtrArray* set) override;
    void   init_table() override;
    void   reduce_cost_fixing(double* pi, int UB, double LB) override;
    void   build_mip() override;
    void   construct_lp_sol_from_rmp(const double*    columns,
                                     const GPtrArray* schedule_sets,
                                     int              num_columns) override;
    void   make_schedule_set_feasible(GPtrArray* set) override;
    void   iterate_zdd() override;
    void   create_dot_zdd(const char* name) override;
    void   print_number_nodes_edges() override;
    void   add_constraint(Job* job, GPtrArray* list, int order) override;
    void   print_num_paths() override;
    void   remove_constraints(int first, int nb_del) override;
    void   update_rows_coeff(int first) override;
    void   insert_constraints_lp(NodeData* pd) override;
    int    get_num_remove_nodes() override;
    int    get_num_remove_edges() override;
    int    get_num_layers() override;
    int    add_constraints() override;
    size_t get_nb_edges() override;
    size_t get_nb_vertices() override;

    inline DdStructure<>* get_decision_diagram() {
        return decision_diagram.get();
    }

    inline int get_nb_removed_edges() { return nb_removed_edges; }

    inline void add_nb_removed_edges() { nb_removed_edges++; }

    inline int* get_take() override { return NULL; };

   private:
    void add_inequality(std::vector<int> v1, std::vector<int> v2);
    void add_inequality(std::vector<int> v1);

    double compute_reduced_cost(const OptimalSolution<>& sol,
                                double*                  pi,
                                double*                  lhs) override;
    double compute_lagrange(const OptimalSolution<>& sol, double* pi) override;

    double compute_subgradient(const OptimalSolution<>& sol,
                               double*                  sub_gradient) override;

    void update_constraints() override {}

    void update_coeff_constraints() override;
};

#endif  // PRICER_SOLVER_BDD_HPP
