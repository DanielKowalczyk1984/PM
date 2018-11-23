#include <wctprivate.h>
#include <PricerSolver.hpp>
#include <iostream>
#include <vector>

template <typename T = double, bool reverse = false>
int construct_sol(wctdata *pd, Optimal_Solution<T> &sol) {
    int          val = 0;
    int          nbset = 1;
    Job *        tmp_j;
    scheduleset *newset = scheduleset_alloc_bis(pd->njobs);
    CCcheck_NULL_3(newset, "Failed to allocate memory newset");

    for (unsigned i = 0; i < sol.jobs->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(sol.jobs, i);
        g_hash_table_add(newset->table, tmp_j);
    }
    newset->jobs = sol.jobs;
    sol.jobs = nullptr;
    newset->e_list = sol.e_list;
    sol.e_list = nullptr;

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

PricerSolver *newSolver(GPtrArray *jobs, GPtrArray *ordered_jobs, int nmachines, int ub, int Hmax) {
    return new PricerSolver(jobs, ordered_jobs, nmachines, ub, Hmax);
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

static double compute_lagrange(Optimal_Solution<double> &sol,
                               double *                  rhs,
                               double *                  pi,
                               int                       nbjobs) {
    double result = 0.0;
    double a = 0.0;
    Job *  tmp_j;

    for (unsigned i = 0; i < sol.jobs->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(sol.jobs, i);
        result += pi[tmp_j->job];
    }

    for (int i = 0; i < nbjobs; ++i) {
        a += rhs[i] * pi[i];
    }

    result = CC_MIN(0.0, sol.cost - result);
    result = rhs[nbjobs] * result;
    result += a;
    return result;
}

static double compute_reduced_cost(Optimal_Solution<double> &sol,
                               double *                  pi,
                               int                       nbjobs) {
    double result = pi[nbjobs];
    Job *  tmp_j;
    int C = 0;

    for (unsigned i = 0; i < sol.jobs->len; ++i) {
        tmp_j = (Job *)g_ptr_array_index(sol.jobs, i);
        C += tmp_j->processingime;
        result += pi[tmp_j->job] - value_Fj(C, tmp_j);
    }

    return sol.obj;
}

int evaluate_nodes(wctdata *pd) {
    int    val = 0;
    int    UB = pd->problem->opt_sol->tw;
    double LB = pd->LP_lower_bound;
    int    nmachines = pd->problem->nmachines;

    pd->solver->evaluate_nodes(pd->pi, UB, LB, nmachines, pd->reduced_cost);

    return val;
}

int calculate_new_ordered_jobs(wctdata *pd){
    int val = 0;

    pd->solver->calculate_new_ordered_jobs();

    return val;
}

int build_solve_mip(wctdata *pd) {
    int val = 0;

    pd->solver->build_mip(pd->x_e);

    return val;
}

void print_number_nodes_edges(wctdata *pd){
    pd->solver->print_number_nodes_edges();
}

void deletePricerSolver(PricerSolver *solver) {
    if (solver) {
        delete solver;
    }
}

int calculate_table(PricerSolver *solver, wctparms *parms) {
    int val = 0;

    solver->init_zdd_table();
    solver->init_table_farkas();

    return val = 0;
}

void iterate_zdd(PricerSolver *solver) { solver->iterate_zdd(); }

void print_number_paths(PricerSolver *solver){ solver->print_number_paths(); }

size_t get_datasize(PricerSolver *solver) { return solver->get_datasize(); }

size_t get_numberrows_zdd(PricerSolver *solver) {
    return solver->get_numberrows_zdd();
}

double get_edge_cost(PricerSolver *solver, int idx) {return solver->get_cost_edge(idx);}

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

int solve_pricing(wctdata *pd, wctparms *parms, int evaluate) {
    int val = 0;


    Optimal_Solution<double> sol = pd->solver->solve_duration_zdd_double(pd->pi);

    // if(pd->iterations%5 == 0 || evaluate) {
    //     pd->reduced_cost = compute_reduced_cost(sol, pd->pi, pd->njobs);
    //     pd->solver->evaluate_nodes(pd->pi, pd->problem->opt_sol->tw, pd->LP_lower_bound_BB , pd->problem->nmachines, pd->reduced_cost);
    // }

    if(sol.obj > 0.000001) {
       val = construct_sol(pd, sol);
       CCcheck_val_2(val, "Failed in construction")
    } else {
        pd->nnewsets = 0;
    }


    if(pd->iterations%pd->njobs == 0) {
        print_number_nodes_edges(pd);
    }

CLEAN:
    return val;
}

int solve_pricing_ti(wctdata *pd, wctparms *parms){
    int val = 0;
    Optimal_Solution<double> sol = pd->solver->dynamic_programming_ti(pd->pi);

    if(sol.obj > 0.000001) {
        val = construct_sol(pd, sol);
        CCcheck_val_2(val, "Failed in construction");
    } else {
        pd->nnewsets = 0;
    }


CLEAN:
    return val;
}


int solve_stab(wctdata *pd, wctparms *parms) {
    int           val = 0;
    PricerSolver *solver = pd->solver;
        double k = 0.0;
        double alpha;
        bool   mispricing = false;
        double result_sep;

    // if((pd->LP_lower_bound_BB - pd->eta_in)/(pd->LP_lower_bound_BB) > 0.5) {
    //     Optimal_Solution<double> sol = solver->solve_weight_zdd_double(pd->pi_out);
    //     result_sep = compute_lagrange(sol, pd->rhs, pd->pi_out, pd->njobs);
    //     pd->reduced_cost = compute_reduced_cost(sol, pd->pi, pd->njobs);
    //     if(pd->reduced_cost >= 0.00001) {
    //         val = construct_sol(pd, sol);
    //         CCcheck_val_2(val, "Failed in construct_sol_stab");
    //         pd->update = 1;
    //     }
    // } else {
        do {
            k += 1.0;
            alpha = CC_MAX(0, 1.0 - k * (1.0 - pd->alpha));
            compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha,
                               pd->pi_in, &(pd->eta_in), pd->pi_out,
                               &(pd->eta_out));
            Optimal_Solution<double> sol;
            if(parms->pricing_solver == bdd_solver) {
                sol = solver->solve_duration_zdd_double(pd->pi_sep);
            } else if (parms->pricing_solver == dp_solver){
                sol = solver->dynamic_programming_ti(pd->pi_sep);
            }
            result_sep = compute_lagrange(sol, pd->rhs, pd->pi_sep, pd->njobs);
            pd->reduced_cost = compute_reduced_cost(sol, pd->pi_out, pd->njobs);

            if (pd->reduced_cost >= 0.00001) {
                val = construct_sol(pd, sol);
                CCcheck_val_2(val, "Failed in construct_sol_stab");
                pd->update = 1;
                mispricing = false;
            } else {
                pd->reduced_cost =compute_reduced_cost(sol, pd->pi_out, pd->njobs);

                if (pd->reduced_cost < 0.000001) {
                    CCcheck_val_2(val, "Failed in construct_sol_stab");
                    mispricing = true;
                    pd->update = 0;
                } else {
                    val = construct_sol(pd, sol);
                    mispricing = false;
                    pd->update = 1;
                }
            }
        } while (mispricing && alpha > 0); /** mispricing check */
    // }


    if (result_sep > pd->eta_in) {
        pd->eta_in = result_sep;
        memcpy(pd->pi_in, pd->pi_sep, sizeof(double) * (pd->njobs + 1));
    }

    if (pd->iterations%pd->njobs == 0) {
        printf(
            "alpha = %f, result of primal bound and Lagragian "
            "bound: out =%f, in = %f\n",
             pd->alpha, pd->eta_out, pd->eta_in);
    }

CLEAN:
    return val;
}

int solve_stab_dynamic(wctdata *pd, wctparms *parms) {
    int           val = 0;
    PricerSolver *solver = pd->solver;
    double        k = 0.0;
    double        alpha;
    double        result_sep;
    double        result_out;
    bool          mispricing = false;

    do {
        k += 1.0;
        alpha = CC_MAX(0.0, 1.0 - k * (1 - pd->alpha));
        compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha,
                           pd->pi_in, &(pd->eta_in), pd->pi_out,
                           &(pd->eta_out));
        Optimal_Solution<double> sol =
            solver->solve_weight_zdd_double(pd->pi_sep);
        result_sep = compute_lagrange(sol, pd->rhs, pd->pi_sep, pd->njobs);

        if (result_sep < pd->eta_sep) {
            val = construct_sol(pd, sol);
            CCcheck_val_2(val, "Failed in construct_sol_stab");
            compute_subgradient(sol, pd->subgradient, pd->rhs, pd->njobs,
                                pd->nmachines);
            adjust_alpha(pd->pi_out, pd->pi_in, pd->subgradient, pd->njobs,
                         alpha);
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

void calculate_edges(PricerSolver *solver, scheduleset *set){
    solver->calculate_edges(set);
}

void g_calculate_edges(gpointer data, gpointer user_data){
    scheduleset *tmp = (scheduleset *) data;
    PricerSolver *solver = (PricerSolver *) user_data;

    solver->calculate_edges(tmp);

}
}
