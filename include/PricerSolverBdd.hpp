#ifndef PRICER_SOLVER_BDD_HPP
#define PRICER_SOLVER_BDD_HPP
// #include <memory>
// #include <range/v3/view/drop.hpp>
// #include <unordered_map>
// #include <vector>
// #include "Instance.h"
// #include "MipGraph.hpp"
// #include "ModelInterface.hpp"
// #include "NodeBdd.hpp"
// #include "NodeBddStructure.hpp"
// #include "OptimalSolution.hpp"
// #include "PricerSolverBase.hpp"

#include <boost/move/utility_core.hpp>           // for move
#include <cstddef>                               // for size_t
#include <memory>                                // for unique_ptr, shared_ptr
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=
#include <range/v3/view/filter.hpp>              // for filter_view
#include <range/v3/view/view.hpp>                // for operator|
#include <utility>                               // for pair
#include <vector>                                // for vector
#include "Instance.h"                            // for Instance
#include "MipGraph.hpp"                          // for MipGraph
#include "ModelInterface.hpp"                    // for OriginalModel
#include "NodeBdd.hpp"                           // for NodeBdd
#include "NodeBddStructure.hpp"                  // for DdStructure
#include "OptimalSolution.hpp"                   // for OptimalSolution
#include "PricerSolverBase.hpp"                  // for PricerSolverBase
struct Interval;
struct Job;
struct NodeData;
struct ScheduleSet;
struct Sol;

class PricerSolverBdd : public PricerSolverBase {
    DdStructure<NodeBdd<double>> decision_diagram;
    size_t                       size_graph;
    size_t                       nb_removed_edges{};
    size_t                       nb_removed_nodes{};

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
    ~PricerSolverBdd() override;

    void                  check_infeasible_arcs();
    void                  topdown_filtering();
    void                  bottum_up_filtering();
    void                  equivalent_paths_filtering();
    [[maybe_unused]] void print_representation_file();
    void                  cleanup_arcs();

    bool refinement_structure(
        const std::vector<std::shared_ptr<ScheduleSet>>& paths) override;
    void enumerate_columns() override;
    void enumerate_columns(double* _pi) override;
    bool evaluate_nodes(double* _pi) override;

    bool check_schedule_set(const std::vector<Job*>& set) override;
    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override = 0;
    virtual double evaluate_rc_arc(NodeBdd<>& n) = 0;
    double         evaluate_rc_low_arc(NodeBdd<>& n);
    virtual void   compute_labels(double* _pi) = 0;

    void remove_layers();
    void remove_edges();
    void remove_layers_init();
    void construct_mipgraph();
    void init_coeff_constraints();
    void init_table();

    // bool check_schedule_set(GPtrArray* set) override;

    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                                    columns,
        const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets)
        override;

    void project_sol_on_original_variables(const Sol& _sol) override;
    std::vector<std::vector<double>>& calculate_job_time() override;
    void   split_job_time(size_t _job, int _time, bool _left) override;
    void   iterate_zdd() override;
    void   create_dot_zdd(const char* name) override;
    void   print_number_nodes_edges() override;
    size_t print_num_paths() override;
    void   remove_constraints(int first, int nb_del) override;
    void   update_rows_coeff(size_t first) override;
    void   insert_constraints_lp(NodeData* pd) override;

    size_t get_num_remove_nodes() override;
    size_t get_num_remove_edges() override;
    int    get_num_layers() override;
    int    add_constraints() override;

    size_t get_nb_edges() override;
    size_t get_nb_vertices() override;
    bool   structure_feasible() override;

    inline DdStructure<NodeBdd<double>>& get_decision_diagram() {
        return decision_diagram;
    }

    inline int get_nb_removed_edges() { return nb_removed_edges; }

    inline void add_nb_removed_edges() { nb_removed_edges++; }

   private:
    void add_inequality(std::vector<int> v1, std::vector<int> v2);
    void add_inequality(std::vector<int> v1);

    double compute_reduced_cost(const OptimalSolution<>& sol,
                                double*                  pi,
                                double*                  lhs) override;
    double compute_lagrange(const OptimalSolution<>&   sol,
                            const std::vector<double>& pi) override;

    double compute_subgradient(const OptimalSolution<>& sol,
                               double*                  sub_gradient) override;

    void update_constraints() override {}

    void update_coeff_constraints() override;
};

#endif  // PRICER_SOLVER_BDD_HPP
