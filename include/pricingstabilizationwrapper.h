#ifndef __PRICINGSTABILIZATIONWRAPPER_H__
#define __PRICINGSTABILIZATIONWRAPPER_H__

// #include "PricingStabilization.hpp"
#include "solver.h"
#include "wctparms.h"

#ifdef __cplusplus
extern "C" {

#endif  // __cplusplus

typedef struct PricingStabilizationBase PricingStabilization;

PricingStabilization* new_pricing_stabilization(PricerSolver* solver,
                                                Parms*        parms);
void   delete_pricing_stabilization(PricingStabilization* pricing_stab_solver);
double call_get_eta_in(PricingStabilization* pricing_stabilization);
double call_get_reduced_cost(PricingStabilization* p);
int    call_stopping_criteria(PricingStabilization* p);
int    call_get_update_stab_center(PricingStabilization* p);
#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // __PRICINGSTABILIZATIONWRAPPER_H__