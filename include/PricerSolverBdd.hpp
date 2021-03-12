#ifndef PRICER_SOLVER_BDD_HPP
#define PRICER_SOLVER_BDD_HPP
#include <memory>
#include <unordered_map>
#include <vector>
#include "Instance.h"
#include "MipGraph.hpp"
#include "ModelInterface.hpp"
#include "NodeBdd.hpp"
#include "NodeBddStructure.hpp"
#include "OptimalSolution.hpp"
#include "PricerSolverBase.hpp"
// #include "wctprivate.h"

class PricerSolverBdd : public PricerSolverBase {
    DdStructure<NodeBdd<double>> decision_diagram;
    size_t                       size_graph;
    int                          nb_removed_edges{};
    int                          nb_removed_nodes{};

    std::vector<std::pair<Job*, Interval*>> ordered_jobs_new;

    MipGraph mip_graph;

    OriginalModel<> original_model;

    // std::unordered_map<int, std::vector<std::weak_ptr<NodeId>>> t_in;
    // std::unordered_map<int, std::vector<std::weak_ptr<NodeId>>> t_out;

    int H_min;
    int H_max;

   public:
    explicit PricerSolverBdd(const Instance& instance);

    PricerSolverBdd(const PricerSolverBdd& src);
    PricerSolverBdd(PricerSolverBdd&&) = default;
    PricerSolverBdd& operator=(PricerSolverBdd&&) = default;
    PricerSolverBdd& operator=(const PricerSolverBdd&) = delete;
    virtual ~PricerSolverBdd() override;
    std::unique_ptr<PricerSolverBase> clone() const override = 0;

    void evaluate_nodes(double* pi, int UB, double LB) override = 0;

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
    void init_table();

    // bool check_schedule_set(GPtrArray* set) override;
    bool check_schedule_set(const std::vector<Job*>& set) override;

    void reduce_cost_fixing(double* pi, int UB, double LB) override;
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                                    columns,
        const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets,
        int                                              num_columns) override;
    // void make_schedule_set_feasible(GPtrArray* set) override;
    void calculate_job_time(std::vector<std::vector<double>>* v) override;
    void split_job_time(int _job, int _time, bool _left) override;
    void iterate_zdd() override;
    void create_dot_zdd(const char* name) override;
    void print_number_nodes_edges() override;
    // void add_constraint(Job* job, GPtrArray* list, int order) override;
    void print_num_paths() override;
    void remove_constraints(int first, int nb_del) override;
    void update_rows_coeff(int first) override;
    void insert_constraints_lp(NodeData* pd) override;

    int get_num_remove_nodes() override;
    int get_num_remove_edges() override;
    int get_num_layers() override;
    int add_constraints() override;

    size_t get_nb_edges() override;
    size_t get_nb_vertices() override;

    inline DdStructure<>& get_decision_diagram() { return decision_diagram; }

    inline int get_nb_removed_edges() { return nb_removed_edges; }

    inline void add_nb_removed_edges() { nb_removed_edges++; }

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
