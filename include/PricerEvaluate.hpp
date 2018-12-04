#include <ForwardBDD.hpp>
#include <ForwardZDD.hpp>

struct ForwardZddCycleDouble : ForwardZddCycle<ForwardZddCycleDouble, double> {
    ForwardZddCycleDouble(double *_pi, int _num_jobs)
    : ForwardZddCycle<ForwardZddCycleDouble, double>(_pi, _num_jobs) {};
    ForwardZddCycleDouble(int _num_jobs)
    : ForwardZddCycle<ForwardZddCycleDouble, double> (_num_jobs) {};
    ForwardZddCycleDouble() : ForwardZddCycle<ForwardZddCycleDouble, double>() {};
};

struct ForwardZddSimpleDouble : ForwardZddSimple<ForwardZddSimpleDouble, double> {
    ForwardZddSimpleDouble(double *_pi, int _num_jobs)
    : ForwardZddSimple<ForwardZddSimpleDouble, double>(_pi, _num_jobs) {};
    ForwardZddSimpleDouble(int _num_jobs)
    : ForwardZddSimple<ForwardZddSimpleDouble, double> (_num_jobs) {};
    ForwardZddSimpleDouble() : ForwardZddSimple<ForwardZddSimpleDouble, double>() {};
};

struct ForwardBddSimpleDouble : ForwardBddSimple<ForwardBddSimpleDouble, double> {
    ForwardBddSimpleDouble(double *_pi, int _num_jobs)
    : ForwardBddSimple<ForwardBddSimpleDouble, double>(_pi, _num_jobs) {};
    ForwardBddSimpleDouble(int _num_jobs)
    : ForwardBddSimple<ForwardBddSimpleDouble, double> (_num_jobs) {};
    ForwardBddSimpleDouble() : ForwardBddSimple<ForwardBddSimpleDouble, double>() {};
};

struct ForwardBddCycleDouble : ForwardBddCycle<ForwardBddCycleDouble, double> {
    ForwardBddCycleDouble(double *_pi, int _num_jobs)
    : ForwardBddCycle<ForwardBddCycleDouble, double>(_pi, _num_jobs) {};
    ForwardBddCycleDouble(int _num_jobs)
    : ForwardBddCycle<ForwardBddCycleDouble, double> (_num_jobs) {};
    ForwardBddCycleDouble() : ForwardBddCycle<ForwardBddCycleDouble, double>() {};
};

