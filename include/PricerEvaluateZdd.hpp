#ifndef PRICER_EVALUATE_ZDD_HPP
#define PRICER_EVALUATE_ZDD_HPP

#include <ForwardZDD.hpp>

struct ForwardZddCycleDouble : ForwardZddCycle<ForwardZddCycleDouble, double> {
    ForwardZddCycleDouble(double* _pi, int _num_jobs)
        : ForwardZddCycle<ForwardZddCycleDouble, double>(_pi, _num_jobs) {};
    explicit ForwardZddCycleDouble(int _num_jobs)
        : ForwardZddCycle<ForwardZddCycleDouble, double> (_num_jobs) {};
    ForwardZddCycleDouble() : ForwardZddCycle<ForwardZddCycleDouble, double>() {};
};

struct ForwardZddSimpleDouble : ForwardZddSimple<ForwardZddSimpleDouble, double> {
    ForwardZddSimpleDouble(double* _pi, int _num_jobs)
        : ForwardZddSimple<ForwardZddSimpleDouble, double>(_pi, _num_jobs) {};
    explicit ForwardZddSimpleDouble(int _num_jobs)
        : ForwardZddSimple<ForwardZddSimpleDouble, double> (_num_jobs) {};
    ForwardZddSimpleDouble() : ForwardZddSimple<ForwardZddSimpleDouble, double>() {};
};

#endif // PRICER_EVALUATE_ZDD_HPP
