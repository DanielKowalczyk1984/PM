#ifndef PRICER_EVALUATE_BDD_HPP
#define PRICER_EVALUATE_BDD_HPP
#include <ForwardBDD.hpp>
#include <BackwardBDD.hpp>
#include <FarkasZDD.hpp>

struct ForwardBddSimpleDouble : ForwardBddSimple<ForwardBddSimpleDouble, double> {
    ForwardBddSimpleDouble() : ForwardBddSimple<ForwardBddSimpleDouble, double>() {};
};

struct ForwardBddCycleDouble : ForwardBddCycle<ForwardBddCycleDouble, double> {
    ForwardBddCycleDouble() : ForwardBddCycle<ForwardBddCycleDouble, double>() {};
};

struct BackwardBddSimpleDouble : BackwardBddSimple<BackwardBddSimpleDouble, double> {
    BackwardBddSimpleDouble() : BackwardBddSimple<BackwardBddSimpleDouble, double>() {};
};

struct BackwardBddFarkasDouble : BackwardBddFarkas<BackwardBddFarkasDouble, double> {
    BackwardBddFarkasDouble() : BackwardBddFarkas<BackwardBddFarkasDouble, double>() {};
};
struct BackwardBddCycleDouble : BackwardBddCycle<BackwardBddCycleDouble, double> {
    BackwardBddCycleDouble() : BackwardBddCycle<BackwardBddCycleDouble, double>() {};
};

#endif // PRICER_EVALUATE_BDD_HPP
