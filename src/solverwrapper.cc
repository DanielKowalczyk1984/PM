#include <wctprivate.h>
#include <PricerSolver.hpp>
#include <iostream>
#include <vector>


template <typename T = double, bool reverse = false>
int construct_sol(wctdata *pd, Optimal_Solution<T> &sol) {
    int               val = 0;
    int               nbset = 1;
    Job *tmp_j;
    scheduleset *     newset = scheduleset_alloc(pd->njobs);
    CCcheck_NULL_3(newset, "Failed to allocate memory newset");

    for(unsigned i = 0; i < sol.jobs->len; ++i) {
        tmp_j = (Job *) g_ptr_array_index(sol.jobs, i);
        g_ptr_array_add(newset->jobs, tmp_j);
        g_hash_table_add(newset->table, tmp_j);
    }


    newset->totwct = sol.cost;
    newset->totweight = sol.C_max;
    pd->newsets = newset;
    pd->nnewsets = 1;
CLEAN:

    if (val) {
        schedulesets_free(&(newset), &(nbset));
    }

    return val;
}

extern "C" {

static double compute_lagrange(Optimal_Solution<double> &sol,
                               double *                  rhs,
                               double *                  pi,
                               int                       nbjobs);

PricerSolver *newSolver(GPtrArray *interval_list, int njobs) {
    return new PricerSolver(interval_list, njobs);
}

PricerSolver *copySolver(PricerSolver *src) { return new PricerSolver(*src); }

void print_dot_file(PricerSolver *solver, char *name) {
    solver->create_dot_zdd(name);
}

void freeSolver(PricerSolver *src) { delete src; }


int solve_farkas_dbl(wctdata *pd) {
    int                      val = 0;
    Optimal_Solution<double> s = pd->solver->solve_farkas_double(pd->pi);

    if (s.obj < -0.00001) {
        // val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->jobarray, s,
        //                     pd->njobs);
        CCcheck_val_2(val, "Failed in constructing jobs");
    } else {
        pd->nnewsets = 0;
    }

CLEAN:
    return val;
}

int solve_weight_dbl_zdd(wctdata *pd) {
    int                      val = 0;
    Optimal_Solution<double> s = pd->solver->solve_weight_zdd_double(pd->pi);

    if (s.obj > 0.00001) {
        val = construct_sol(pd, s);
        CCcheck_val_2(val, "Failed in construction")
    } else {
        pd->nnewsets = 0;
    }

CLEAN:
    return val;
}

void deletePricerSolver(PricerSolver *solver) {
    if (solver) {
        delete solver;
    }
}

int calculate_table(PricerSolver *solver, wctparms *parms) {
    int val = 0;

    solver->init_zdd_table();

    switch (parms->construct) {
        case yes_construct:
            break;

        case no_construct:
            solver->init_table_farkas();
            break;
    }

    return val = 0;
}

int add_conflict_constraints(PricerSolver *solver,
                             wctparms *    parms,
                             int *         elist_same,
                             int           ecount_same,
                             int *         elist_differ,
                             int           ecount_differ) {
    int val = 0;


    solver->init_zdd_conflict_solver(elist_same, ecount_same,
                                             elist_differ, ecount_differ);

    return val;
}

void iterate_zdd(PricerSolver *solver) { solver->iterate_zdd(); }

int free_conflict_constraints(PricerSolver *solver,
                              wctparms *    parms,
                              int           ecount_same,
                              int           ecount_differ) {
    int val = 0;

    solver->free_zdd_solver(ecount_same, ecount_differ);

    return val;
}

size_t get_datasize(PricerSolver *solver) { return solver->zdd->size(); }

size_t get_numberrows_zdd(PricerSolver *solver) {
    return solver->zdd->root().row();
}

int add_one_conflict(
    PricerSolver *solver, wctparms *parms, Job *v1, Job *v2, int same) {
    int val = 0;

    solver->init_zdd_one_conflict(v1->job, v2->job, same);

    return val;
}

int init_tables(PricerSolver *solver) {
    int val = 0;
    solver->init_tables();
    return val;
}

static void compute_subgradient(Optimal_Solution<double> &sol,
                                double *                  sub_gradient,
                                double *                  rhs,
                                int                       nbjobs,
                                int                       nbmachines) {
    fill_dbl(sub_gradient, nbjobs, 1.0);
    sub_gradient[nbjobs] = nbmachines;

    // for (auto &v : sol.jobs) {
    //     sub_gradient[v] -= (double)nbmachines * 1.0;
    // }
}

static void adjust_alpha(double *pi_out,
                         double *pi_in,
                         double *subgradient,
                         int     nbjobs,
                         double &alpha) {
    double sum = 0.0;

    for (int i = 0; i <= nbjobs; ++i) {
        sum += subgradient[i] * (pi_out[i] - pi_in[i]);
    }

    if (sum > 0) {
        alpha = alpha + (1 - alpha) * 0.05;
    } else {
        alpha = CC_MAX(0, alpha - 0.05);
    }
}

static void compute_pi_eta_sep(int     vcount,
                               double *pi_sep,
                               double *eta_sep,
                               double  alpha,
                               double *pi_in,
                               double *eta_in,
                               double *pi_out,
                               double *eta_out) {
    int    i;
    double beta = 1.0 - alpha;

    for (i = 0; i <= vcount; ++i) {
        pi_sep[i] = alpha * pi_in[i] + beta * pi_out[i];
    }

    *eta_sep = alpha * (*eta_in) + beta * (*eta_out);
}

int solve_pricing(wctdata *pd, wctparms *parms) {
    int val = 0;

    val = solve_weight_dbl_zdd(pd);
    CCcheck_val_2(val, "Failed solve_weight_dbl_zdd");

CLEAN:
    return val;
}

static double compute_lagrange(Optimal_Solution<double> &sol,
                               double *                  rhs,
                               double *                  pi,
                               int                       nbjobs) {
    double            result = 0.0;
    double a = 0.0;
    Job *tmp_j;


    for(unsigned i = 0; i < sol.jobs->len; ++i) {
        tmp_j = (Job *) g_ptr_array_index(sol.jobs, i);
        result += pi[tmp_j->job];
    }

    for(int i = 0; i < nbjobs; ++i) {
        a += rhs[i]*pi[i];
    }

    result = CC_MIN(0.0, sol.cost - result);
    result = rhs[nbjobs] * result;
    result += a;
    return result;
}

int solve_stab(wctdata *pd, wctparms *parms) {
    int                      val = 0;
    int                      heading_in = 0;
    PricerSolver *           solver = pd->solver;
    double result;
    heading_in =
        (pd->eta_in == 0.0)
            ? 1
            : !(CC_OURABS((pd->eta_out -
                     pd->eta_in) / (pd->eta_in)) < 4.0);

    if (heading_in) {
        Optimal_Solution<double> sol =  solver->solve_weight_zdd_double(pd->pi);
        result = compute_lagrange(sol, pd->rhs, pd->pi, pd->njobs);
        if (result > pd->eta_in) {
            pd->eta_in = result;
            memcpy(pd->pi_in, pd->pi, sizeof(double) * (pd->njobs + 1));
        }

        val = construct_sol(pd, sol);
        CCcheck_val_2(val, "Failed in construct solution");
    } else {
        double k = 0.0;
        double alpha;
        bool   mispricing = false;
        double result_sep;
        double result_out;

        do {
            k += 1.0;
            alpha = CC_MAX(0, 1.0 - k * (1.0 - pd->alpha));
            compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha,
                               pd->pi_in, &(pd->eta_in), pd->pi_out,
                               &(pd->eta_out));
            Optimal_Solution<double> sol = solver->solve_weight_zdd_double(pd->pi_sep);
            result_sep = compute_lagrange(sol, pd->rhs, pd->pi_sep, pd->njobs);

            if (result_sep < pd->eta_sep) {
                val = construct_sol(pd, sol);
                CCcheck_val_2(val, "Failed in construct_sol_stab");
                pd->update = 1;
                mispricing = false;
            } else {
                result_out = compute_lagrange(sol, pd->rhs, pd->pi_out, pd->njobs);

                if (result_out < pd->eta_out) {
                    val = construct_sol(pd, sol);
                    CCcheck_val_2(val, "Failed in construct_sol_stab");
                    mispricing = false;
                    pd->update = 0;
                } else {
                    mispricing = true;
                }
            }
        } while (mispricing && alpha > 0); /** mispricing check */

        if (result_sep > pd->eta_in) {
            pd->eta_in = result_sep;
            memcpy(pd->pi_in, pd->pi_sep, sizeof(double) * (pd->njobs + 1));
        }
    }

    if (dbg_lvl() > 0) {
        printf(
            "heading = %d, alpha = %f, result of primal bound and Lagragian "
            "bound: out =%f, in = %f\n",
            heading_in, pd->alpha, pd->eta_out, pd->eta_in);
    }

CLEAN:
    return val;
}

int solve_stab_dynamic(wctdata *pd, wctparms *parms) {
    int                      val = 0;
    PricerSolver *           solver = pd->solver;
    double                   k = 0.0;
    double                   alpha;
    double                   result_sep;
    double                   result_out;
    bool                     mispricing = false;

    do {
        k += 1.0;
        alpha = CC_MAX(0.0, 1.0 - k * (1 - pd->alpha));
        compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha,pd->pi_in, &(pd->eta_in), pd->pi_out,&(pd->eta_out));
        Optimal_Solution<double> sol = solver->solve_weight_zdd_double(pd->pi_sep);
        result_sep = compute_lagrange(sol, pd->rhs, pd->pi_sep, pd->njobs);

        if (result_sep < pd->eta_sep) {
            val = construct_sol(pd, sol);
            CCcheck_val_2(val, "Failed in construct_sol_stab");
            compute_subgradient(sol, pd->subgradient, pd->rhs, pd->njobs, pd->nmachines);
            adjust_alpha(pd->pi_out, pd->pi_in, pd->subgradient, pd->njobs, alpha);
            pd->alpha = alpha;
            pd->update = 1;
            mispricing = false;
        } else {
            result_out = compute_lagrange(sol, pd->rhs, pd->pi_out, pd->njobs);

            if (result_out < pd->eta_out) {
                val = construct_sol(pd, sol);
                CCcheck_val_2(val, "Failed in construct_sol_stab");
                mispricing = false;
                pd->update = 0;
            } else {
                mispricing = true;
            }
        }
    } while (mispricing && alpha > 0.0);

    if (result_sep > pd->eta_in) {
        pd->eta_in = result_sep;
        memcpy(pd->pi_in, pd->pi_sep, sizeof(double) * (pd->njobs + 1));
    }

    if (dbg_lvl() > 0) {
        printf(
            " alpha = %f, result of primal bound and Lagragian bound: out =%f, "
            "in = %f\n",
            pd->alpha, pd->eta_out, pd->eta_in);
    }

CLEAN:
    return val;
}
}
