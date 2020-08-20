#ifndef PRICER_EVALUATE_BDD_HPP
#define PRICER_EVALUATE_BDD_HPP
#include <BackwardBDD.hpp>
#include <FarkasZDD.hpp>
#include <ForwardBDD.hpp>

using ForwardBddSimpleDouble = ForwardBddSimple<>;
using ForwardBddCycleDouble = ForwardBddCycle<>;
using BackwardBddSimpleDouble = BackwardBddSimple<>;
using BackwardBddFarkasDouble = BackwardBddFarkas<>;
using BackwardBddCycleDouble = BackwardBddCycle<>;

#endif  // PRICER_EVALUATE_BDD_HPP
