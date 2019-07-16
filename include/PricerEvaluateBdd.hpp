#ifndef PRICER_EVALUATE_BDD_HPP
#define PRICER_EVALUATE_BDD_HPP
#include <ForwardBDD.hpp>
#include <BackwardBDD.hpp>

struct ForwardBddSimpleDouble : ForwardBddSimple<ForwardBddSimpleDouble, double> {
    ForwardBddSimpleDouble(double *_pi, int _num_jobs)
    : ForwardBddSimple<ForwardBddSimpleDouble, double>(_pi, _num_jobs) {};
    explicit ForwardBddSimpleDouble(int _num_jobs)
    : ForwardBddSimple<ForwardBddSimpleDouble, double> (_num_jobs) {};
    ForwardBddSimpleDouble() : ForwardBddSimple<ForwardBddSimpleDouble, double>() {};
};

struct ForwardBddCycleDouble : ForwardBddCycle<ForwardBddCycleDouble, double> {
    ForwardBddCycleDouble(double *_pi, int _num_jobs)
    : ForwardBddCycle<ForwardBddCycleDouble, double>(_pi, _num_jobs) {};
    explicit ForwardBddCycleDouble(int _num_jobs)
    : ForwardBddCycle<ForwardBddCycleDouble, double> (_num_jobs) {};
    ForwardBddCycleDouble() : ForwardBddCycle<ForwardBddCycleDouble, double>() {};
};

struct BackwardBddSimpleDouble : BackwardBddSimple<BackwardBddSimpleDouble, double> {
    BackwardBddSimpleDouble(double *_pi, int _num_jobs)
    : BackwardBddSimple<BackwardBddSimpleDouble, double>(_pi, _num_jobs) {};
    explicit BackwardBddSimpleDouble(int _num_jobs)
    : BackwardBddSimple<BackwardBddSimpleDouble, double> (_num_jobs) {};
    BackwardBddSimpleDouble() : BackwardBddSimple<BackwardBddSimpleDouble, double>() {};
};

struct BackwardBddCycleDouble : BackwardBddCycle<BackwardBddCycleDouble, double> {
    BackwardBddCycleDouble(double *_pi, int _num_jobs)
    : BackwardBddCycle<BackwardBddCycleDouble, double>(_pi, _num_jobs) {};
    explicit BackwardBddCycleDouble(int _num_jobs)
    : BackwardBddCycle<BackwardBddCycleDouble, double> (_num_jobs) {};
    BackwardBddCycleDouble() : BackwardBddCycle<BackwardBddCycleDouble, double>() {};
};

#endif // PRICER_EVALUATE_BDD_HPP
