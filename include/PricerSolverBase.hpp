#ifndef PRICER_SOLVER_BASE_HPP
#define PRICER_SOLVER_BASE_HPP

#include <gurobi_c++.h>
#include <solution.h>
#include <MIP_defs.hpp>
#include <OptimalSolution.hpp>
#include <cstddef>
#include <memory>
#include "ModelInterface.hpp"
#include "scheduleset.h"
#include "wctprivate.h"

struct PricerSolverBase {
   public:
    GPtrArray*                jobs;
    int                       nb_jobs;
    int                       num_machines;
    GPtrArray*                ordered_jobs;
    int                       nb_layers;
    std::string               problem_name;
    std::unique_ptr<GRBEnv>   env;
    std::unique_ptr<GRBModel> model;
    ReformulationModel        reformulation_model;
    bool                      is_integer_solution;
    /**
     * Default constructors
     */
    PricerSolverBase(GPtrArray* _jobs, int _num_machines, const char* p_name);
    PricerSolverBase(GPtrArray*  _jobs,
                     int         _num_machines,
                     GPtrArray*  _ordered_jobs,
                     const char* p_name);
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

    /**
     * init_table
     */
    virtual void init_table() = 0;

    /**
     * Pricing Algorithm
     */
    virtual OptimalSolution<double> pricing_algorithm(double* _pi) = 0;
    virtual OptimalSolution<double> farkas_pricing(double* _pi) = 0;

    /**
     * Reduced cost fixing
     */

    virtual void reduce_cost_fixing(double* pi, int UB, double LB) = 0;
    virtual void evaluate_nodes(double* pi, int UB, double LB) = 0;

    /** Original Mip formulation */
    virtual void build_mip() = 0;
    virtual void construct_lp_sol_from_rmp(const double*    columns,
                                           const GPtrArray* schedule_sets,
                                           int              num_columns) = 0;
    virtual void represent_solution(Solution* sol) = 0;
    virtual void project_solution(Solution* sol) = 0;

    /**
     * Constraint on the solver
     */

    virtual void add_constraint(Job* job, GPtrArray* list, int order) = 0;
    virtual void insert_constraints_lp(NodeData* pd) = 0;
    virtual void add_constraints(){

    };

    virtual void update_coeff_constraints() = 0;

    /**
     * Some getters
     */
    virtual void iterate_zdd() = 0;
    virtual void print_num_paths() = 0;

    virtual int    get_num_remove_nodes() = 0;
    virtual int    get_num_remove_edges() = 0;
    virtual int    get_num_layers() = 0;
    virtual size_t get_nb_vertices() = 0;
    virtual size_t get_nb_edges() = 0;
    virtual bool   check_schedule_set(GPtrArray* set) = 0;
    virtual void   make_schedule_set_feasible(GPtrArray* set) = 0;
    /**
     * Some printing functions
     */
    virtual void   create_dot_zdd(const char* name) = 0;
    virtual void   print_number_nodes_edges() = 0;
    virtual int*   get_take() = 0;
    virtual int    get_int_attr_model(enum MIP_Attr);
    virtual double get_dbl_attr_model(enum MIP_Attr);

    /**
     * Constraints auxilary functions
     */

    virtual void   update_constraints() = 0;
    virtual double compute_reduced_cost(const OptimalSolution<>& sol,
                                        double*                  pi,
                                        double*                  lhs);
    virtual double compute_lagrange(const OptimalSolution<>& sol, double* pi);

    virtual double compute_subgradient(const OptimalSolution<>& sol,
                                       double*                  pi,
                                       double*                  subgradient){};

    inline void set_is_integer_solution(bool _is_solution) {
        is_integer_solution = _is_solution;
    }

    inline bool get_is_integer_solution() { return is_integer_solution; }

    inline ReformulationModel* get_reformulation_model() {
        return &reformulation_model;
    }
    // virtual void compute_lhs_coeff(GArray *lhs_coeff, ScheduleSet* set) = 0;
};

#endif  // PRICER_SOLVER_BASE_HPP
