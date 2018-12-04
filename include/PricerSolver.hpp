#ifndef INCLUDE_PRICERSOLVER_HPP
#define INCLUDE_PRICERSOLVER_HPP

#include <PricerEvaluate.hpp>
#include <PricerConstruct.hpp>
#include <tdzdd/DdStructure.hpp>
#include <scheduleset.h>
#include <DurationZDD.hpp>
#include <DurationBDD.hpp>
using namespace std;

struct PricerSolverBase {
protected:
    GPtrArray *jobs;
    int njobs;
    GPtrArray *ordered_jobs;
    int nlayers;

    tdzdd::DdStructure<2> *zdd;
    tdzdd::DdStructure<2> *dd;

    size_t nb_nodes_bdd;
    size_t nb_nodes_zdd;

    int nb_removed_edges;
    int nb_removed_nodes;
public:
    /**
     * Default constructors
     */
    PricerSolverBase(GPtrArray *_jobs);
    PricerSolverBase(GPtrArray *_jobs, GPtrArray *_ordered_jobs);
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
     * InitTable
     */
    virtual void InitTable() = 0;

    /**
     * Pricing Algorithm
     */
     virtual Optimal_Solution<double> pricing_algorithm(double *_pi) = 0;

     /**
      * Reduced cost fixing
      */
     void calculate_new_ordered_jobs();
     void calculate_edges(scheduleset *set);
     void evaluate_nodes(double *pi, int UB, double LB, int nmachines,
                         double reduced_cost);

     /**
      * Some getters
      */
     void IterateZdd();
     void PrintNumberPaths();

     int get_remove();
     size_t get_datasize();
     size_t get_numberrows_zdd();
     double get_cost_edge(int idx);

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
    PricerSolverBdd(GPtrArray *_jobs, GPtrArray *_ordered_jobs);
    void InitTable();
};

class PricerSolverBddSimple : public PricerSolverBdd {
private:
    ForwardBddSimpleDouble evaluator;
public:
    PricerSolverBddSimple(GPtrArray *_jobs, GPtrArray *_ordered_jobs);
    Optimal_Solution<double> pricing_algorithm(double *_pi);
};

class PricerSolverBddCycle : public PricerSolverBdd {
private:
    ForwardBddCycleDouble evaluator;
public:
    PricerSolverBddCycle(GPtrArray *_jobs, GPtrArray *_ordered_jobs);
    Optimal_Solution<double> pricing_algorithm(double *_pi);
};


class PricerSolverZdd : public PricerSolverBase {
protected:
    tdzdd::DataTable<ForwardZddNode<double>> table;
public:
    PricerSolverZdd(GPtrArray *_jobs, GPtrArray *ordered_jobs);
    void InitTable();
};

class PricerSolverZddSimple : public PricerSolverZdd {
private:
    ForwardZddSimpleDouble evaluator;
public:
    PricerSolverZddSimple(GPtrArray *_jobs, GPtrArray *_ordered_jobs);

    /**
     * Pricing Algorithm
     */
    Optimal_Solution<double> pricing_algorithm(double * _pi);
};

class PricerSolverCycle : public PricerSolverZdd {
private:
    ForwardZddCycleDouble evaluator;
public:
    PricerSolverCycle(GPtrArray *_jobs, GPtrArray *_ordered_jobs);

    /**
     * Pricing Algorithm
     */
    Optimal_Solution<double> pricing_algorithm(double * _pi);

};

class PricerSolverSimpleDp : public PricerSolverBase {
private:
    int Hmax;
public:
    PricerSolverSimpleDp(GPtrArray *_jobs, int _Hmax);
    void InitTable();

    Optimal_Solution<double> pricing_algorithm(double *_pi);
};

#endif  // INCLUDE_PRICERSOLVER_HPP
