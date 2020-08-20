#ifndef PRICER_EVALUATE_ZDD_HPP
#define PRICER_EVALUATE_ZDD_HPP

#include "BackwardZDD.hpp"
#include "ForwardZDD.hpp"

// struct ForwardZddCycleDouble : ForwardZddCycle<ForwardZddCycleDouble, double>
// {
//     ForwardZddCycleDouble(double* _pi, int _num_jobs)
//         : ForwardZddCycle<ForwardZddCycleDouble, double>(_pi, _num_jobs){};
//     explicit ForwardZddCycleDouble(int _num_jobs)
//         : ForwardZddCycle<ForwardZddCycleDouble, double>(_num_jobs){};
//     ForwardZddCycleDouble()
//         : ForwardZddCycle<ForwardZddCycleDouble, double>(){};
// };

using ForwardZddCycleDouble = ForwardZddCycle<>;

// struct ForwardZddSimpleDouble
//     : ForwardZddSimple<ForwardZddSimpleDouble, double> {
//     ForwardZddSimpleDouble(double* _pi, int _num_jobs)
//         : ForwardZddSimple<ForwardZddSimpleDouble, double>(_pi, _num_jobs){};
//     explicit ForwardZddSimpleDouble(int _num_jobs)
//         : ForwardZddSimple<ForwardZddSimpleDouble, double>(_num_jobs){};
//     ForwardZddSimpleDouble()
//         : ForwardZddSimple<ForwardZddSimpleDouble, double>(){};
// };

using ForwardZddSimpleDouble = ForwardZddSimple<>;

// struct BackwardZddSimpleDouble
//     : BackwardZddSimple<BackwardZddSimpleDouble, double> {
//     BackwardZddSimpleDouble(double* _pi, int _num_jobs)
//         : BackwardZddSimple<BackwardZddSimpleDouble, double>(_pi,
//         _num_jobs){};
//     explicit BackwardZddSimpleDouble(int _num_jobs)
//         : BackwardZddSimple<BackwardZddSimpleDouble, double>(_num_jobs){};
//     BackwardZddSimpleDouble()
//         : BackwardZddSimple<BackwardZddSimpleDouble, double>(){};
// };

using BackwardZddSimpleDouble = BackwardZddSimple<>;

// struct BackwardZddCycleDouble
//     : BackwardZddCycle<BackwardZddCycleDouble, double> {
//     BackwardZddCycleDouble(double* _pi, int _num_jobs)
//         : BackwardZddCycle<BackwardZddCycleDouble, double>(_pi, _num_jobs){};
//     explicit BackwardZddCycleDouble(int _num_jobs)
//         : BackwardZddCycle<BackwardZddCycleDouble, double>(_num_jobs){};
//     BackwardZddCycleDouble()
//         : BackwardZddCycle<BackwardZddCycleDouble, double>(){};
// };

using BackwardZddCycleDouble = BackwardZddCycle<>;

#endif  // PRICER_EVALUATE_ZDD_HPP
