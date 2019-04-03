#ifndef INCLUDE_PRICERSOLVER_HPP
#define INCLUDE_PRICERSOLVER_HPP

#include <tdzdd/DdStructure.hpp>
#include <NodeBddStructure.hpp>
#include <PricerEvaluate.hpp>
#include <scheduleset.h>
#include <MipGraph.hpp>

struct PricerSolverBase {
protected:
    GPtrArray *jobs;
    int njobs;
    int num_machines;
    GPtrArray *ordered_jobs;
    int nlayers;

    std::unique_ptr<DdStructure<>> decision_diagram;

    size_t size_graph;

    int nb_removed_edges;
    int nb_removed_nodes;

    MipGraph g;
    std::unique_ptr<GRBEnv> env;
    std::unique_ptr<GRBModel> model;
public:
    /**
     * Default constructors
     */
    PricerSolverBase(GPtrArray *_jobs, int _num_machines);
    PricerSolverBase(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    /**
     * Copy constructor
     */
    PricerSolverBase(const PricerSolverBase &other);

    /**
     * Move Constructor
     */
    PricerSolverBase(PricerSolverBase &&other) noexcept;

    /**
     * Move Constructor
     */
    PricerSolverBase &operator=(const PricerSolverBase &other);

    /**
     * Move assignment operator
     */
    PricerSolverBase &operator=(PricerSolverBase &&other) noexcept;

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
     virtual Optimal_Solution<double> pricing_algorithm(double *_pi) = 0;

     /**
      * Reduced cost fixing
      */
     void calculate_new_ordered_jobs(double *pi, int UB, double LB);
     void remove_layers();
     void remove_edges();
     void calculate_edges(scheduleset *set);
     virtual void evaluate_nodes(double *pi, int UB, double LB);

     /** Construct MipGraph */
     void construct_mipgraph();
     void build_mip();
     void disjunctive_inequality(double *x, solution *sol);
     void represent_solution(solution *sol);
     double* project_solution(solution *sol);

     /**
      * Some getters
      */
     void iterate_zdd();
     void print_num_paths();

     int get_num_remove_nodes();
     int get_num_remove_edges();
     int get_num_layers();
     size_t get_size_graph();
     size_t get_datasize();
     double get_cost_edge(int idx);
     void construct_lp_sol_from_rmp(const double *columns, const GPtrArray* schedule_sets, int num_columns, double *x);
     bool check_schedule_set(scheduleset *set);

     /**
      * Some printing functions
      */
     void create_dot_zdd(const char *name);
     void print_number_nodes_edges();
};

class PricerSolverBdd : public PricerSolverBase {
protected:
    tdzdd::DataTable<Node<double>> table;
public:
    PricerSolverBdd(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    void init_table() override;
    virtual void evaluate_nodes(double *pi, int UB, double LB) override = 0;
};

class PricerSolverBddSimple : public PricerSolverBdd {
private:
    ForwardBddSimpleDouble evaluator;
    BackwardBddSimpleDouble reversed_evaluator;
public:
    PricerSolverBddSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    Optimal_Solution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverBddCycle : public PricerSolverBdd {
private:
    ForwardBddCycleDouble evaluator;
    BackwardBddCycleDouble reversed_evaluator;
public:
    PricerSolverBddCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    Optimal_Solution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;        
};

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble reversed_evaluator;

public:
    PricerSolverBddBackwardSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);
    Optimal_Solution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
private:
    BackwardBddCycleDouble evaluator;
    ForwardBddCycleDouble reversed_evaluator;
public:
    PricerSolverBddBackwardCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);

    Optimal_Solution<double> pricing_algorithm(double *_pi) override;
    void compute_labels(double *_pi);
    void evaluate_nodes(double *pi, int UB, double LB) override ;    
};

class PricerSolverSimpleDp : public PricerSolverBase {
private:
    int Hmax;
    std::unique_ptr<Job*[]> A;
    std::unique_ptr<double[]> F;
    
public:
    PricerSolverSimpleDp(GPtrArray *_jobs, int _num_machines, int _Hmax);
    ~PricerSolverSimpleDp();
    void init_table() override;

    Optimal_Solution<double> pricing_algorithm(double *_pi) override;
};

class PricerSolverArcTimeDp : public PricerSolverBase {
private:
    int Hmax;
    int n;
    boost::unordered_set <Job*> **graph;
    std::vector<Job*> vector_jobs;
    Job j0;
    double **F;
    Job ***A;
    int **B;
    int **p_matrix;

    typedef boost::unordered_set<Job*>::iterator job_iterator;
public:
    PricerSolverArcTimeDp(GPtrArray *_jobs, int _num_machines, int _Hmax);
    ~PricerSolverArcTimeDp();
    void init_table() override;

    Optimal_Solution<double> pricing_algorithm(double *_pi) override;
    int delta1(const int &i,const int &j, const int &t) {
        Job *tmp_i = vector_jobs[i];
        Job *tmp_j = vector_jobs[j];
        return (value_Fj(t, tmp_i) + value_Fj(t + tmp_j->processing_time, tmp_j))
            - (value_Fj(t + tmp_j->processing_time - tmp_i->processing_time, tmp_j)
                + value_Fj(t + tmp_j->processing_time, tmp_i) );
    }

    void remove_arc(const int &i, const int &j, const int &t) {
        Job *tmp_i = vector_jobs[i];
        // auto it = graph[j][t].find(tmp_i);
        graph[j][t].erase(tmp_i) ;
    }

    int delta2(const int &j, const int &t) {
        Job *tmp_j = vector_jobs[j];
        return value_Fj(t, tmp_j) - value_Fj(t + 1, tmp_j);
    }
};


// class PricerSolverZdd : public PricerSolverBase {
// protected:
//     tdzdd::DataTable<ForwardZddNode<double>> table;
//     ForwardZddNode<double>& child(tdzdd::NodeId const & id);
// public:
//     PricerSolverZdd(GPtrArray *_jobs, int _num_machines, GPtrArray *ordered_jobs);
//     void init_table() override;
// };

// class PricerSolverZddSimple : public PricerSolverZdd {
// private:
//     ForwardZddSimpleDouble evaluator;
// public:
//     PricerSolverZddSimple(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);

//     /**
//      * Pricing Algorithm
//      */
//     Optimal_Solution<double> pricing_algorithm(double * _pi) override;
// };

// class PricerSolverCycle : public PricerSolverZdd {
// private:
//     ForwardZddCycleDouble evaluator;
// public:
//     PricerSolverCycle(GPtrArray *_jobs, int _num_machines, GPtrArray *_ordered_jobs);

//     /**
//      * Pricing Algorithm
//      */
//     Optimal_Solution<double> pricing_algorithm(double * _pi) override;

// };

#endif  // INCLUDE_PRICERSOLVER_HPP
