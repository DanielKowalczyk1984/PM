#ifndef PRICER_SOLVER_ZDD_HPP
#define PRICER_SOLVER_ZDD_HPP
#include <NodeBddStructure.hpp>                  // for DdStructure
#include <cstddef>                               // for size_t
#include <memory>                                // for unique_ptr, shared_ptr
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=
#include <range/v3/view/filter.hpp>              // for filter_view
#include <range/v3/view/view.hpp>                // for operator|
#include <utility>                               // for pair
#include <vector>                                // for vector
#include "Instance.h"                            // for Instance
#include "MipGraph.hpp"                          // for MipGraph
#include "OptimalSolution.hpp"                   // for OptimalSolution
#include "PricerSolverBase.hpp"                  // for PricerSolverBase
#include "ZddNode.hpp"                           // for NodeZdd
struct Interval;                                 // lines 36-36
struct Job;                                      // lines 37-37
struct NodeData;                                 // lines 38-38
struct Column;                                   // lines 39-39
class PricerSolverZdd : public PricerSolverBase {
   public:
    std::unique_ptr<DdStructure<NodeZdd<double>>> decision_diagram;

    size_t size_graph{};
    size_t nb_removed_edges{};
    size_t nb_removed_nodes{};

    // GPtrArray*                              ordered_jobs;
    std::vector<std::pair<Job*, Interval*>> ordered_jobs_new;

    MipGraph            mip_graph;
    std::vector<double> lp_x;
    std::vector<double> solution_x;

    explicit PricerSolverZdd(const Instance& instance);
    PricerSolverZdd(PricerSolverZdd&&) = default;
    PricerSolverZdd& operator=(PricerSolverZdd&&) = default;
    PricerSolverZdd& operator=(const PricerSolverZdd&) = delete;

    PricerSolverZdd(const PricerSolverZdd& src)
        : PricerSolverBase(src),
          // decision_diagram(new
          // DdStructure<NodeZdd<>>(*src.decision_diagram)),
          size_graph(src.size_graph),
          nb_removed_edges(src.nb_removed_edges),
          nb_removed_nodes(src.nb_removed_nodes),
          //   ordered_jobs(src.ordered_jobs),
          ordered_jobs_new(src.ordered_jobs_new),
          mip_graph(src.mip_graph) {}

    ~PricerSolverZdd() override = default;
    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return nullptr;
    };

    void init_table();

    void remove_layers();
    void remove_edges();

    void construct_mipgraph();
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                               columns,
        const std::vector<std::shared_ptr<Column>>& schedule_sets) override;
    bool    check_schedule_set(const std::vector<Job*>& set) override;
    void    iterate_zdd() override;
    void    create_dot_zdd(const char* name) override;
    size_t  get_nb_edges() override;
    size_t  get_nb_vertices() override;
    cpp_int print_num_paths() override;
    void    remove_layers_init();

    OptimalSolution<double> farkas_pricing(double* pi) override;

    // void add_constraint(Job* job, GPtrArray* list, int order) override;

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

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
