#ifndef PRICER_SOLVER_BASE_HPP
#define PRICER_SOLVER_BASE_HPP
#include <gurobi_c++.h>         // for GRBModel
#include <cstddef>              // for size_t
#include <functional>           // for reference_wrapper, ref
#include <memory>               // for shared_ptr, unique_ptr
#include <string>               // for string
#include <utility>              // for ref
#include <vector>               // for vector
#include "Instance.h"           // for Instance
#include "MIP_defs.hpp"         // for MIP_Attr
#include "ModelInterface.hpp"   // for BddCoeff, ReformulationModel
#include "OptimalSolution.hpp"  // for OptimalSolution
#include "Solution.hpp"         // for Sol
struct Job;
struct NodeData;     // lines 15-15
struct ScheduleSet;  // lines 16-16

struct PricerSolverBase {
   public:
    const std::vector<std::shared_ptr<Job>>& jobs;

    size_t convex_constr_id;
    size_t convex_rhs;

    std::string problem_name;

    std::shared_ptr<GRBEnv> env;
    GRBModel                model;

    ReformulationModel reformulation_model;

    bool is_integer_solution;

    double constLB;
    double UB;

    std::vector<std::vector<double>> x_bar;
    std::vector<std::vector<double>> z_bar;

    static const std::shared_ptr<GRBEnv> genv;
    static constexpr double              EPS_SOLVER = 1e-6;
    static constexpr double              RC_FIXING = 1e-4;
    static constexpr int                 ALIGN = 60;
    static constexpr int                 ALIGN_HALF = 60;

    std::vector<BddCoeff> lp_sol;
    /**
     * Default constructors
     */
    explicit PricerSolverBase(const Instance& instance);

    /**
     * Copy constructor
     */
    PricerSolverBase(const PricerSolverBase& other);

    /**
     * Move Constructor
     */
    PricerSolverBase(PricerSolverBase&& other) noexcept;

    /**
     * Move Constructor
     */
    PricerSolverBase& operator=(const PricerSolverBase& other);

    /**
     * Move assignment operator
     */
    PricerSolverBase& operator=(PricerSolverBase&& other) noexcept;

    /**
     * Destructor
     */
    virtual ~PricerSolverBase();

    [[nodiscard]] virtual std::unique_ptr<PricerSolverBase> clone() const = 0;

    /**
     * Pricing Algorithm
     */
    virtual OptimalSolution<double> pricing_algorithm(double* _pi) = 0;
    virtual OptimalSolution<double> farkas_pricing(double* _pi) = 0;

    /**
     * Reduced cost fixing
     */

    virtual bool evaluate_nodes(double* pi) = 0;
    virtual bool refinement_structure(
        [[maybe_unused]] const std::vector<std::shared_ptr<ScheduleSet>>&
            paths) {
        return false;
    };
    virtual void enumerate_columns(){};
    virtual void enumerate_columns([[maybe_unused]] double* _pi){};

    /** Original Mip formulation */
    virtual void build_mip() = 0;
    bool         evaluate_mip_model();
    virtual void construct_lp_sol_from_rmp(
        const double*                                    columns,
        const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets) = 0;

    bool compute_sub_optimal_duals(
        const double*                                    columns,
        const std::vector<std::shared_ptr<ScheduleSet>>& schedule_sets);
    virtual void project_sol_on_original_variables(const Sol& _sol) {
        _sol.print_solution();
    };

    /**
     * Constraint on the solver
     */

    virtual void   insert_constraints_lp(NodeData* pd) = 0;
    virtual int    add_constraints();
    virtual void   remove_constraints(int first, int nb_del);
    virtual void   update_rows_coeff(size_t first);
    virtual size_t get_size_data() { return 0UL; };

    virtual void                              update_coeff_constraints() = 0;
    virtual std::vector<std::vector<double>>& calculate_job_time() {
        return x_bar;
    };

    std::pair<std::vector<std::vector<double>>&,
              std::vector<std::vector<double>>&>
    get_pair_x() {
        return std::make_pair(std::ref(x_bar), std::ref(z_bar));
    }

    virtual void split_job_time([[maybe_unused]] size_t _job,
                                [[maybe_unused]] int    _time,
                                [[maybe_unused]] bool   _left) {}

    /**
     * Some getters
     */
    virtual void   iterate_zdd() = 0;
    virtual size_t print_num_paths() = 0;
    double         get_UB();
    void           update_UB(double _ub);

    virtual size_t get_num_remove_nodes() = 0;
    virtual size_t get_num_remove_edges() = 0;
    virtual int    get_num_layers() = 0;
    virtual size_t get_nb_vertices() = 0;
    virtual size_t get_nb_edges() = 0;
    virtual bool   structure_feasible() { return true; }
    virtual bool   check_schedule_set(
          [[maybe_unused]] const std::vector<Job*>& set) {
        return true;
    };
    virtual void make_schedule_set_feasible(
        [[maybe_unused]] std::vector<Job*>& set){};
    /**
     * Some printing functions
     */
    virtual void   create_dot_zdd(const char* name) = 0;
    virtual void   print_number_nodes_edges() = 0;
    virtual int    get_int_attr_model(enum MIP_Attr);
    virtual double get_dbl_attr_model(enum MIP_Attr);

    /**
     * Constraints auxilary functions
     */

    virtual void   update_constraints() = 0;
    virtual double compute_reduced_cost(const OptimalSolution<>& sol,
                                        double*                  pi,
                                        double*                  lhs);

    virtual double compute_reduced_cost_simple(const OptimalSolution<>& sol,
                                               double*                  pi);

    virtual double compute_lagrange(const OptimalSolution<>&   sol,
                                    const std::vector<double>& pi);

    virtual double compute_subgradient(const OptimalSolution<>& sol,
                                       double*                  subgradient);

    inline void set_is_integer_solution(bool _is_solution) {
        is_integer_solution = _is_solution;
    }

    inline bool get_is_integer_solution() { return is_integer_solution; }

    inline ReformulationModel* get_reformulation_model() {
        return &reformulation_model;
    }

    inline std::vector<BddCoeff>& get_lp_sol();

    void calculate_constLB(double* pi);
};

#endif  // PRICER_SOLVER_BASE_HPP
