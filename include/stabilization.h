#ifndef STABILIZATION_H
#define STABILIZATION_H

#include <wctprivate.h>
#include <OptimalSolution.hpp>

#ifdef __cplusplus
extern "C" {
#endif

void   update_alpha(NodeData* pd);
void   update_alpha_misprice(NodeData* pd);
int    is_stabilized(NodeData* pd);
int    calculate_dualdiffnorm(NodeData* pd);
int    calculate_beta(NodeData* pd);
int    calculate_hybridfactor(NodeData* pd);
int    update_hybrid(NodeData* pd);
int    update_node(NodeData* pd);
int    get_dual_row(NodeData* pd, int i);
int    row_getDual(NodeData* pd, int i);
int    update_subgradientproduct(NodeData* pd);
int    update_stabcenter(const OptimalSolution<double>& sol, NodeData* pd);
double compute_dual(NodeData* pd, int i);
double compute_lagrange(OptimalSolution<double>& sol, double* rhs, double* pi,
                        int nbjobs);

#ifdef __cplusplus
}
#endif

#endif  // STABILIZATION_H
