#ifndef INCLUDE_PRICERSOLVER_HPP
#define INCLUDE_PRICERSOLVER_HPP

#include <PricerEvaluate.hpp>
#include <PricerConstruct.hpp>
#include <tdzdd/DdStructure.hpp>
#include <scheduleset.h>
#include <DurationZDD.hpp>
using namespace std;

struct PricerSolver {
  private:
    int nmachines;
    int ub;
    int Hmax;
    int njobs;
    int nlayers;
    GPtrArray *jobs;
    GPtrArray *ordered_jobs;
    tdzdd::DdStructure<2>                    *zdd;
    tdzdd::DdStructure<2>                    *dd;
    tdzdd::DataTable<PricerWeightZDD<double>> zdd_table;
    tdzdd::DataTable<PricerDurationZDD<double>> zdd_duration_table;
    tdzdd::DataTable<PricerFarkasZDD<double>> farkas_table;
    tdzdd::DataTable<PricerInfoBDD<double> > dd_table;
    std::vector<shared_ptr<edge<double>>> edges;
    WeightZDDdouble evaluator_weight;
    DurationZDDdouble evaluator;
    int nb_removed_edges;
    int nb_removed_nodes;
    size_t nb_nodes_bdd;
    size_t nb_nodes_zdd;
    GRBEnv *env;
    GRBModel *model;


  public:

    /** Default Constructor */
    PricerSolver(GPtrArray *_jobs, GPtrArray *_ordered_jobs, int _nmachines,
                 int _ub, int _Hmax);

    /** Copy constructor */
    PricerSolver(const PricerSolver &other);

    /** Move Constructor */
    PricerSolver(PricerSolver &&other) noexcept;

    /** Move Constructor */
    PricerSolver &operator=(const PricerSolver &other);

    /** Move assignment operator */
    PricerSolver &operator=(PricerSolver &&other) noexcept;

    /** Destructor */
    ~PricerSolver() noexcept;

    int get_remove();
    size_t get_datasize();
    size_t get_numberrows_zdd();
    double get_cost_edge(int idx);

    
    /**
     * Reduce the number of nodes and edges in the decision diagram
     */
    void calculate_new_ordered_jobs();
    void calculate_edges(scheduleset *set);
    void evaluate_nodes(double *pi, int UB, double LB, int nmachines,
                        double reduced_cost);

    /**
     * Some printing functions
     */
    void create_dot_zdd(const char *name);
    void print_number_nodes_edges();

    /**
     * Init tables 
     */
    void init_zdd_table();
    void init_bdd_table();
    void init_table_farkas();
    void init_zdd_duration_table();
    void init_tables();

    /**
     * Iteration of all paths
     */
    void iterate_zdd();
    void print_number_paths();

    /**
     * Pricing Algorithms
     */
    Optimal_Solution<double> dynamic_programming_ti(double *pi);
    Optimal_Solution<double> solve_duration_bdd_double(double *pi);
    Optimal_Solution<double> solve_weight_zdd_double(double *pi);
    Optimal_Solution<double> solve_farkas_double(double *pi);
    Optimal_Solution<double> solve_duration_zdd_double(double *pi);
    
    /**
     * build the reduced extented formulation   
     */
    void build_mip(double *x_e);
};

#endif  // INCLUDE_PRICERSOLVER_HPP
