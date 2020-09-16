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
void delete_pricing_stabilization(PricingStabilization* pricing_stabilization);
double call_get_eta_in(PricingStabilization* pricing_stabilization);
double call_get_reduced_cost(PricingStabilization* pricing_stabilization);
int    call_stopping_criteria(PricingStabilization* pricing_stabilization);
int    call_get_update_stab_center(PricingStabilization* pricing_stabilization);
void   call_update_duals(PricingStabilization* pricing_stabilization);
void   call_remove_constraints(PricingStabilization* pricing_stabilization,
                               int                   first,
                               int                   nb_del);
#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // __PRICINGSTABILIZATIONWRAPPER_H__