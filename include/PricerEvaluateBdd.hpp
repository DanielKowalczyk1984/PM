#ifndef PRICER_EVALUATE_BDD_HPP
#define PRICER_EVALUATE_BDD_HPP
#include <BackwardBDD.hpp>
#include <FarkasZDD.hpp>
#include <ForwardBDD.hpp>
#include "ModelInterface.hpp"

// struct ForwardBddSimpleDouble
//     : ForwardBddSimple<ForwardBddSimpleDouble, double> {
//     ForwardBddSimpleDouble()
//         : ForwardBddSimple<ForwardBddSimpleDouble, double>(){};
//     ForwardBddSimpleDouble(OriginalModel<>* model)
//         : ForwardBddSimple<ForwardBddSimpleDouble, double>(model){};
// };

using ForwardBddSimpleDouble = ForwardBddSimple<>;

// struct ForwardBddCycleDouble : ForwardBddCycle<ForwardBddCycleDouble, double>
// {
//     ForwardBddCycleDouble()
//         : ForwardBddCycle<ForwardBddCycleDouble, double>(){};
//     ForwardBddCycleDouble(OriginalModel<>* model)
//         : ForwardBddCycle<ForwardBddCycleDouble, double>(model){};
// };

using ForwardBddCycleDouble = ForwardBddCycle<>;

// struct BackwardBddSimpleDouble
//     : BackwardBddSimple<BackwardBddSimpleDouble, double> {
//     BackwardBddSimpleDouble()
//         : BackwardBddSimple<BackwardBddSimpleDouble, double>(){};
//     BackwardBddSimpleDouble(OriginalModel<>* model)
//         : BackwardBddSimple<BackwardBddSimpleDouble, double>(model){};
// };

using BackwardBddSimpleDouble = BackwardBddSimple<>;

// struct BackwardBddFarkasDouble
//     : BackwardBddFarkas<BackwardBddFarkasDouble, double> {
//     BackwardBddFarkasDouble()
//         : BackwardBddFarkas<BackwardBddFarkasDouble, double>(){};
// };

using BackwardBddFarkasDouble = BackwardBddFarkas<>;

// struct BackwardBddCycleDouble
//     : BackwardBddCycle<BackwardBddCycleDouble, double> {
//     BackwardBddCycleDouble()
//         : BackwardBddCycle<BackwardBddCycleDouble, double>(){};
//     BackwardBddCycleDouble(OriginalModel<>* model)
//         : BackwardBddCycle<BackwardBddCycleDouble, double>(model){};
// };

using BackwardBddCycleDouble = BackwardBddCycle<>;

#endif  // PRICER_EVALUATE_BDD_HPP
