#ifndef STABILIZATION_H
#define STABILIZATION_H

#include <wctprivate.h>
#include <OptimalSolution.hpp>

#ifdef __cplusplus
extern "C" {
#endif

void update_alpha(wctdata *pd);
void update_alpha_misprice(wctdata *pd);
int is_stabilized(wctdata *pd);
int calculate_dualdiffnorm(wctdata *pd);
int calculate_beta(wctdata *pd);
int calculate_hybridfactor(wctdata *pd);
int update_hybrid(wctdata *pd);
int update_node(wctdata *pd);
int get_dual_row(wctdata *pd, int i);
int row_getDual(wctdata *pd, int i);
int update_subgradientproduct(wctdata *pd);
int update_stabcenter(const Optimal_Solution<double> & sol, wctdata *pd);
double compute_dual(wctdata *pd, int i);
double compute_lagrange(Optimal_Solution<double> &sol,
                               double *                  rhs,
                               double *                  pi,
                               int                       nbjobs);


#ifdef __cplusplus
}
#endif

#endif // STABILIZATION_H
