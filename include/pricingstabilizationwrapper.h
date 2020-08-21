#ifndef __PRICINGSTABILIZATIONWRAPPER_H__
#define __PRICINGSTABILIZATIONWRAPPER_H__

// #include "PricingStabilization.hpp"
#include "solver.h"

#ifdef __cplusplus
extern "C" {

#endif  // __cplusplus

typedef struct PricingStabilization PricingStabilization;

PricingStabilization* new_pricing_stabilization(PricerSolver* solver);
void   delete_pricing_stabilization(PricingStabilization* pricing_stab_solver);
double get_diff_eta_out_eta_in(PricingStabilization* pricing_stabilization);
double call_C_get_diff_eta(PricingStabilization* p);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // __PRICINGSTABILIZATIONWRAPPER_H__