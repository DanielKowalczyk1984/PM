#include <stabilization.h>
#include <wctprivate.h>
#include <cstdio>
#include "PricerSolverBase.hpp"
#include "lp.h"
#include "scheduleset.h"

template <typename T = double>
int construct_sol(NodeData* pd, OptimalSolution<T>* sol) {
    int          val = 0;
    int          nb_set = 1;
    ScheduleSet* newset = scheduleset_alloc_bis(pd->nb_jobs);
    CCcheck_NULL_3(newset, "Failed to allocate memory newset");

    for (unsigned i = 0; i < sol->jobs->len; ++i) {
        Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol->jobs, i));
        g_hash_table_add(newset->table, tmp_j);
    }
    newset->job_list = sol->jobs;
    sol->jobs = nullptr;
    newset->edge_list = nullptr;

    newset->total_weighted_completion_time = sol->cost;
    newset->total_processing_time = sol->C_max;
    pd->newsets = newset;
    pd->nb_new_sets = 1;

CLEAN:
    if (val) {
        schedulesets_free(&(newset), &(nb_set));
    }

    return val;
}

extern "C" {

double compute_lagrange(OptimalSolution<double>& sol,
                        double*                  rhs,
                        double*                  pi,
                        int                      nb_jobs) {
    double result = 0.0;
    double a = 0.0;

    for (unsigned i = 0; i < sol.jobs->len; ++i) {
        Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol.jobs, i));
        result += pi[tmp_j->job];
    }

    for (int i = 0; i < nb_jobs; ++i) {
        a += rhs[i] * pi[i];
    }

    result = CC_MIN(0.0, sol.cost - result);
    result = -rhs[nb_jobs] * result;
    result += a;

    return result;
}

// static double compute_reduced_cost(OptimalSolution<double>& sol, double* pi,
//                                    int nb_jobs) {
//     double result = sol.cost + pi[nb_jobs];

//     for (guint i = 0; i < sol.jobs->len; ++i) {
//         Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol.jobs, i));
//         result -= pi[tmp_j->job];
//     }

//     return result;
// }

void update_alpha(NodeData* pd) {
    if (pd->subgradientproduct > 0.0) {
        pd->alpha = CC_MAX(0, pd->alpha - 0.1);
    } else {
        pd->alpha = CC_MIN(0.99, pd->alpha + (1 - pd->alpha) * 0.1);
    }
}

void update_alpha_misprice(NodeData* pd) {
    pd->k++;
    pd->alphabar = CC_MAX(0.0, 1 - pd->k * (1.0 - pd->alpha));
}

int is_stabilized(NodeData* pd) {
    if (pd->in_mispricing_schedule) {
        return pd->alphabar > 0.0;
    }
    return pd->alpha > 0.0;
}

int calculate_dualdiffnorm(NodeData* pd) {
    int     val = 0;
    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* pi_in = &g_array_index(pd->pi_in, double, 0);

    pd->dualdiffnorm = 0.0;

    for (int i = 0; i <= pd->nb_jobs; ++i) {
        double dualdiff = SQR(pi_in[i] - pi_out[i]);
        if (dualdiff > 0.00001) {
            pd->dualdiffnorm += dualdiff;
        }
    }

    pd->dualdiffnorm = SQRT(pd->dualdiffnorm);

    return val;
}

int calculate_beta(NodeData* pd) {
    int     val = 0;
    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* pi_in = &g_array_index(pd->pi_in, double, 0);
    double* subgradient_in = &g_array_index(pd->subgradient_in, double, 0);

    pd->beta = 0.0;
    for (int i = 0; i < pd->nb_rows; ++i) {
        double dualdiff = ABS(pi_out[i] - pi_in[i]);
        double product = dualdiff * ABS(subgradient_in[i]);

        if (product > 0.000001) {
            pd->beta += product;
        }
    }

    if (pd->subgradientnorm > 0.00001) {
        pd->beta = pd->beta / (pd->subgradientnorm * pd->dualdiffnorm);
    }

    return val;
}

int calculate_hybridfactor(NodeData* pd) {
    int val = 0;

    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* pi_in = &g_array_index(pd->pi_in, double, 0);
    double* subgradient_in = &g_array_index(pd->subgradient_in, double, 0);

    double aux_norm = 0.0;
    for (int i = 0; i < pd->nb_rows; ++i) {
        double aux_double =
            SQR((pd->beta - 1.0) * (pi_out[i] - pi_in[i]) +
                pd->beta * (subgradient_in[i] * pd->dualdiffnorm /
                            pd->subgradientnorm));
        if (aux_double > 0.00001) {
            aux_norm += aux_double;
        }
    }
    aux_norm = SQRT(aux_norm);

    pd->hybridfactor = ((1 - pd->alpha) * pd->dualdiffnorm) / aux_norm;

    return val;
}

int update_hybrid(NodeData* pd) {
    int val = 0;

    if (pd->hasstabcenter && !pd->in_mispricing_schedule && pd->alpha > 0.0) {
        calculate_dualdiffnorm(pd);
        calculate_beta(pd);
        calculate_hybridfactor(pd);
    }

    return val;
}

int update_node(NodeData* pd) {
    int val = 0;
    if (pd->node_stab != pd->id) {
        pd->node_stab = pd->id;
        pd->k = 0;
        pd->alpha = 0.0;
        pd->hasstabcenter = 0;
        pd->eta_in = 0.0;
        pd->in_mispricing_schedule = 0;
    }

    return val;
}

double compute_dual(NodeData* pd, int i) {
    double usedalpha = pd->alpha;
    double usedbeta = pd->beta;

    double* pi_in = &g_array_index(pd->pi_in, double, 0);
    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* subgradient_in = &g_array_index(pd->subgradient_in, double, 0);

    if (pd->in_mispricing_schedule) {
        usedalpha = pd->alphabar;
        usedbeta = 0.0;
    }

    if (pd->hasstabcenter && (usedbeta == 0.0 || usedalpha == 0.0)) {
        return usedalpha * pi_in[i] + (1.0 - usedalpha) * pi_out[i];
    } else if (pd->hasstabcenter && usedbeta > 0.0) {
        double dual =
            pi_in[i] +
            pd->hybridfactor *
                (pd->beta * (pi_in[i] + subgradient_in[i] * pd->dualdiffnorm /
                                            pd->subgradientnorm) +
                 (1.0 - pd->beta) * pi_out[i] - pi_in[i]);
        return CC_MAX(dual, 0.0);
    }

    return pi_out[i];
}

int row_getDual(NodeData* pd, int i) {
    int val = 0;
    assert(i < pd->nb_rows);

    g_array_index(pd->pi_sep, double, i) = compute_dual(pd, i);

    return val;
}

static void compute_subgradient(const OptimalSolution<double>& sol,
                                NodeData*                      pd) {
    double* subgradient_in = &g_array_index(pd->subgradient_in, double, 0);
    double* rhs = &g_array_index(pd->rhs, double, 0);
    fill_dbl(subgradient_in, pd->nb_rows, 1.0);
    subgradient_in[pd->nb_jobs] = 0.0;

    for (guint i = 0; i < sol.jobs->len; i++) {
        Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol.jobs, i));
        subgradient_in[tmp_j->job] += rhs[pd->nb_jobs] * 1.0;
    }

    pd->subgradientnorm = 0.0;

    for (int i = 0; i < pd->nb_rows; ++i) {
        double sqr = SQR(subgradient_in[i]);

        if (sqr > 0.00001) {
            pd->subgradientnorm += sqr;
        }
    }

    pd->subgradientnorm = SQRT(pd->subgradientnorm);
}

int update_subgradientproduct(NodeData* pd) {
    int val = 0;

    double* pi_in = &g_array_index(pd->pi_in, double, 0);
    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* subgradient_in = &g_array_index(pd->subgradient_in, double, 0);

    pd->subgradientproduct = 0.0;
    for (int i = 0; i < pd->nb_rows; ++i) {
        pd->subgradientproduct -= (pi_out[i] - pi_in[i]) * subgradient_in[i];
    }
    // printf("subgradientproduct %f\n", pd->subgradientproduct);

    return val;
}

int update_stabcenter(const OptimalSolution<double>& sol, NodeData* pd) {
    int val = 0;

    double* pi_sep = &g_array_index(pd->pi_sep, double, 0);
    double* pi_in = &g_array_index(pd->pi_in, double, 0);

    if (pd->eta_sep > pd->eta_in) {
        memcpy(pi_in, pi_sep, (pd->nb_rows) * sizeof(double));
        compute_subgradient(sol, pd);
        pd->eta_in = pd->eta_sep;
        pd->hasstabcenter = 1;
    }

    return val;
}

static void adjust_alpha(double* pi_out,
                         double* pi_in,
                         double* subgradient,
                         int     nb_jobs,
                         double& alpha) {
    double sum = 0.0;

    for (int i = 0; i < nb_jobs; ++i) {
        sum += subgradient[i] * (pi_out[i] - pi_in[i]);
    }

    if (sum > 0) {
        alpha = CC_MAX(0, alpha - 0.1);
    } else {
        alpha = CC_MIN(0.9, alpha + (1 - alpha) * 0.05);
    }
}

static void compute_pi_eta_sep(int     nb_constr,
                               double* pi_sep,
                               double* eta_sep,
                               double  alpha,
                               double* pi_in,
                               double* eta_in,
                               double* pi_out,
                               double* eta_out) {
    int    i;
    double beta = 1.0 - alpha;

    for (i = 0; i < nb_constr; ++i) {
        pi_sep[i] = alpha * pi_in[i] + beta * pi_out[i];
    }

    *eta_sep = alpha * (*eta_in) + beta * (*eta_out);
}

int solve_pricing(NodeData* pd) {
    int val = 0;

    OptimalSolution<double> sol;
    pd->update = 0;
    double* pi = &g_array_index(pd->pi, double, 0);
    double* lhs = &g_array_index(pd->lhs_coeff, double, 0);

    sol = pd->solver->pricing_algorithm(pi);
    pd->reduced_cost = pd->solver->compute_reduced_cost(sol, pi, lhs);

    if (pd->reduced_cost < -1e-6) {
        val = construct_sol(pd, &sol);
        pd->update = 1;
        CCcheck_val_2(val, "Failed in construction")
    } else {
        pd->nb_new_sets = 0;
    }

CLEAN:
    return val;
}

int solve_stab(NodeData* pd) {
    int           val = 0;
    PricerSolver* solver = pd->solver;
    double        k = 0.0;
    double        alpha;
    bool          mispricing = true;
    double        result_sep;
    double*       pi_sep = &g_array_index(pd->pi_sep, double, 0);
    double*       pi_out = &g_array_index(pd->pi_out, double, 0);
    double*       pi_in = &g_array_index(pd->pi_in, double, 0);
    double*       lhs_coeff = &g_array_index(pd->lhs_coeff, double, 0);
    pd->update = 0;

    do {
        k += 1.0;
        alpha =
            pd->hasstabcenter ? CC_MAX(0, 1.0 - k * (1.0 - pd->alpha)) : 0.0;
        compute_pi_eta_sep(pd->nb_rows, pi_sep, &(pd->eta_sep), alpha, pi_in,
                           &(pd->eta_in), pi_out, &(pd->eta_out));
        OptimalSolution<double> sol;
        sol = solver->pricing_algorithm(pi_sep);

        result_sep = pd->solver->compute_lagrange(sol, pi_sep);
        pd->reduced_cost =
            pd->solver->compute_reduced_cost(sol, pi_out, lhs_coeff);

        if (pd->reduced_cost < -1e-6) {
            val = construct_sol(pd, &sol);
            CCcheck_val_2(val, "Failed in construct_sol_stab");
            pd->update = 1;
            mispricing = false;
        }
    } while (mispricing && alpha > 0); /** mispricing check */

    if (result_sep > pd->eta_in) {
        pd->hasstabcenter = 1;
        pd->eta_in = result_sep;
        memcpy(pi_in, pi_sep, sizeof(double) * (pd->nb_rows));
        pd->update_stab_center = 1;
    } else {
        pd->update_stab_center = 0;
    }

    if (pd->iterations % pd->nb_jobs == 0) {
        printf(
            "alpha = %f, result of primal bound and Lagragian "
            "bound: out =%f, in = %f\n",
            pd->alpha, pd->eta_out, pd->eta_in);
    }

CLEAN:
    return val;
}

int solve_stab_dynamic(NodeData* pd) {
    int           val = 0;
    PricerSolver* solver = pd->solver;
    double        k = 0.0;
    double        alpha;
    double        result_sep;
    bool          mispricing = true;
    pd->update = 0;
    double* pi_sep = &g_array_index(pd->pi_sep, double, 0);
    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* pi_in = &g_array_index(pd->pi_in, double, 0);
    double* subgradient_in = &g_array_index(pd->subgradient_in, double, 0);
    double* lhs_coeff = &g_array_index(pd->lhs_coeff, double, 0);

    do {
        k += 1.0;
        alpha =
            pd->hasstabcenter ? CC_MAX(0.0, 1.0 - k * (1 - pd->alpha)) : 0.0;
        compute_pi_eta_sep(pd->nb_rows, pi_sep, &(pd->eta_sep), alpha, pi_in,
                           &(pd->eta_in), pi_out, &(pd->eta_out));
        OptimalSolution<double> sol;
        sol = solver->pricing_algorithm(pi_sep);
        result_sep = pd->solver->compute_lagrange(sol, pi_sep);
        pd->reduced_cost =
            pd->solver->compute_reduced_cost(sol, pi_out, lhs_coeff);

        if (pd->reduced_cost <= -0.000001) {
            compute_subgradient(sol, pd);
            adjust_alpha(pi_out, pi_in, subgradient_in, pd->nb_jobs, alpha);
            val = construct_sol(pd, &sol);
            CCcheck_val_2(val, "Failed in construct_sol_stab");
            pd->alpha = alpha;
            pd->update = 1;
            mispricing = false;
        }
    } while (mispricing && alpha > 0.0);

    if (result_sep > pd->eta_in) {
        pd->hasstabcenter = 1;
        pd->eta_in = result_sep;
        memcpy(pi_in, pi_sep, sizeof(double) * (pd->nb_rows));
    }

    if (pd->iterations % pd->nb_jobs == 0) {
        printf(
            " alpha = %f, result of primal bound and Lagragian bound: out =%f, "
            "in = %f\n",
            pd->alpha, pd->eta_out, pd->eta_in);
    }

CLEAN:
    return val;
}

int solve_stab_hybrid(NodeData* pd) {
    int           val = 0;
    PricerSolver* solver = pd->solver;
    pd->update = 0;
    bool    stabilized = false;
    double* pi_sep = &g_array_index(pd->pi_sep, double, 0);
    double* pi_out = &g_array_index(pd->pi_out, double, 0);
    double* lhs_coeff = &g_array_index(pd->lhs_coeff, double, 0);

    do {
        update_node(pd);
        update_hybrid(pd);

        stabilized = is_stabilized(pd);

        for (int i = 0; i < pd->nb_rows; ++i) {
            pi_sep[i] = compute_dual(pd, i);
        }

        OptimalSolution<double> sol;
        sol = solver->pricing_algorithm(pi_sep);

        pd->eta_sep = pd->solver->compute_lagrange(sol, pi_sep);
        pd->reduced_cost =
            pd->solver->compute_reduced_cost(sol, pi_out, lhs_coeff);

        update_stabcenter(sol, pd);

        if (pd->reduced_cost < 0.00001) {
            if (pd->in_mispricing_schedule) {
                pd->in_mispricing_schedule = 0;
            }
            update_subgradientproduct(pd);
            update_alpha(pd);
            construct_sol(pd, &sol);
            pd->update = 1;
        } else {
            if (stabilized) {
                pd->in_mispricing_schedule = 1;
                update_alpha_misprice(pd);
            } else {
                pd->in_mispricing_schedule = 0;
            }
        }
    } while (pd->in_mispricing_schedule && stabilized);

    if (pd->iterations % pd->nb_jobs == 0) {
        printf(
            " alpha = %f, result of primal bound and Lagragian bound: out =%f, "
            "in = %f\n",
            pd->alpha, pd->eta_out, pd->eta_in);
    }

    return val;
}

int solve_farkas_dbl(NodeData* pd) {
    int                     val = 0;
    OptimalSolution<double> s =
        pd->solver->farkas_pricing(&g_array_index(pd->pi, double, 0));
    pd->update = 0;

    if (s.obj < 1e-6) {
        val = construct_sol(pd, &s);
        wctlp_write(pd->RMP, "RMP.lp");

        CCcheck_val_2(val, "Failed in constructing jobs");
        pd->update = 1;
    } else {
        pd->nb_new_sets = 0;
    }

CLEAN:
    return val;
}
}
