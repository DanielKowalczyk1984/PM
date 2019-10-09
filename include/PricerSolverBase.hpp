#ifndef PRICER_SOLVER_BASE_HPP
#define PRICER_SOLVER_BASE_HPP

#include <OptimalSolution.hpp>
#include <solution.h>

struct PricerSolverBase {
   protected:
    GPtrArray*  jobs;
    int         nb_jobs;
    int         num_machines;
    GPtrArray*  ordered_jobs;
    int         nb_layers;
    std::string problem_name;

   public:
    /**
     * Default constructors
     */
    PricerSolverBase(GPtrArray* _jobs, int _num_machines, const char* p_name);
    PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                     GPtrArray* _ordered_jobs, const char* p_name);
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

    virtual void disjunctive_inequality(double* x, Solution* sol) = 0;

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

    /**
     * Some printing functions
     */
    virtual void create_dot_zdd(const char* name) = 0;
    virtual void print_number_nodes_edges() = 0;
};

#endif  // PRICER_SOLVER_BASE_HPP
