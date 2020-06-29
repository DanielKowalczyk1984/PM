#ifndef PRICER_EVALUATE_BDD_HPP
#define PRICER_EVALUATE_BDD_HPP
#include <ForwardBDD.hpp>
#include <BackwardBDD.hpp>
#include <FarkasZDD.hpp>
#include "ModelInterface.hpp"

struct ForwardBddSimpleDouble : ForwardBddSimple<ForwardBddSimpleDouble, double> {
    ForwardBddSimpleDouble() : ForwardBddSimple<ForwardBddSimpleDouble, double>() {};
    ForwardBddSimpleDouble(OriginalModel<>*model) : ForwardBddSimple<ForwardBddSimpleDouble, double>(model) {};
};

struct ForwardBddCycleDouble : ForwardBddCycle<ForwardBddCycleDouble, double> {
    ForwardBddCycleDouble() : ForwardBddCycle<ForwardBddCycleDouble, double>() {};
    ForwardBddCycleDouble(OriginalModel<>* model) : ForwardBddCycle<ForwardBddCycleDouble, double>(model) {};
};

struct BackwardBddSimpleDouble : BackwardBddSimple<BackwardBddSimpleDouble, double> {
    BackwardBddSimpleDouble() : BackwardBddSimple<BackwardBddSimpleDouble, double>() {};
    BackwardBddSimpleDouble(OriginalModel<>* model) : BackwardBddSimple<BackwardBddSimpleDouble, double>(model) {};
};

struct BackwardBddFarkasDouble : BackwardBddFarkas<BackwardBddFarkasDouble, double> {
    BackwardBddFarkasDouble() : BackwardBddFarkas<BackwardBddFarkasDouble, double>() {};
};
struct BackwardBddCycleDouble : BackwardBddCycle<BackwardBddCycleDouble, double> {
    BackwardBddCycleDouble() : BackwardBddCycle<BackwardBddCycleDouble, double>() {};
    BackwardBddCycleDouble(OriginalModel<>* model) : BackwardBddCycle<BackwardBddCycleDouble, double>(model) {};
};

#endif // PRICER_EVALUATE_BDD_HPP
