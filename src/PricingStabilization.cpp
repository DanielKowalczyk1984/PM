#include "PricingStabilization.hpp"
#include <fmt/core.h>
#include <cmath>
#include "wctparms.h"

/**
 * @brief Construct a new Pricing Stabilization Base:: Pricing Stabilization
 * Base object
 *
 * @param _solver
 */
PricingStabilizationBase::PricingStabilizationBase(PricerSolverBase* _solver)
    : solver(_solver) {}

/**
@brief solve the pricing problem without stabilization
 *
 * @param _eta_out
 * @param _pi_out
 * @param _lhs
 */
void PricingStabilizationBase::solve(double  _eta_out,
                                     double* _pi_out,
                                     double* _lhs) {
    solver->calculate_constLB(_pi_out);
    sol = std::move(solver->pricing_algorithm(_pi_out));
    reduced_cost = solver->compute_reduced_cost(sol, _pi_out, _lhs);
    iterations++;
    if (iterations % 10 == 0 && _eta_out - 1e-2 <= solver->UB) {
        solver->calculate_constLB(_pi_out);
        solver->evaluate_nodes(_pi_out);
    }
}

OptimalSolution<>& PricingStabilizationBase::get_sol() {
    return sol;
}

double PricingStabilizationBase::get_reduced_cost() {
    return reduced_cost;
}

bool PricingStabilizationBase::get_update_stab_center() {
    return update_stab_center;
}

double PricingStabilizationBase::get_eta_in() {
    return 0.0;
}

int PricingStabilizationBase::stopping_criteria() {
    return reduced_cost < -1e-6;
}

PricingStabilizationBase::~PricingStabilizationBase() {}

/**
 * @brief Wentgnes stabilization technique
 *
 */

/**
 * @brief Construct a new Pricing Stabilization Stat:: Pricing Stabilization
 * Stat object
 *
 * @param _solver
 */
PricingStabilizationStat::PricingStabilizationStat(PricerSolverBase* _solver)
    : PricingStabilizationBase(_solver),
      pi_in(_solver->convex_constr_id + 1, 0.0),
      pi_out(_solver->convex_constr_id + 1),
      pi_sep(_solver->convex_constr_id + 1) {}

PricingStabilizationStat::~PricingStabilizationStat() {}

/**
@brief Solve pricing problem with Wentgnes stabilization
 *
 * @param _eta_out
 * @param _pi_out
 * @param _lhs_coeff
 */
void PricingStabilizationStat::solve(double  _eta_out,
                                     double* _pi_out,
                                     double* _lhs_coeff) {
    k = 0.0;
    bool   mispricing = true;
    double result_sep;
    update = 0;

    std::copy(_pi_out,
              _pi_out + solver->reformulation_model.get_nb_constraints(),
              pi_out.begin());
    eta_out = _eta_out;
    iterations++;
    update_stab_center = false;

    do {
        k += 1.0;
        alphabar = hasstabcenter ? CC_MAX(0, 1.0 - k * (1.0 - alpha)) : 0.0;
        compute_pi_eta_sep(alphabar);
        OptimalSolution<> aux_sol =
            std::move(solver->pricing_algorithm(pi_sep.data()));

        result_sep = solver->compute_lagrange(aux_sol, pi_sep.data());
        reduced_cost =
            solver->compute_reduced_cost(aux_sol, pi_out.data(), _lhs_coeff);

        if (reduced_cost < -1e-6) {
            sol = std::move(aux_sol);
            update = 1;
            mispricing = false;
        }
    } while (mispricing && alphabar > 0); /** mispricing check */

    if (result_sep > eta_in) {
        hasstabcenter = 1;
        eta_in = result_sep;
        pi_in = pi_sep;
        update_stab_center = true;
        if (iterations % 10 == 0 && eta_out - 1e-2 <= solver->UB) {
            solver->calculate_constLB(pi_sep.data());
            solver->evaluate_nodes(pi_sep.data());
        }
    }

    if (iterations % solver->convex_constr_id == 0) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }
}

double PricingStabilizationStat::get_eta_in() {
    return eta_in;
}

int PricingStabilizationStat::stopping_criteria() {
    return (std::abs(eta_out - eta_in) >= 1e-4);
}

/**
 * Dynamic stabilization technique
 */

PricingStabilizationDynamic::PricingStabilizationDynamic(
    PricerSolverBase* _solver)
    : PricingStabilizationStat(_solver),
      subgradient(_solver->convex_constr_id + 1, 0.0)

{}

PricingStabilizationDynamic::~PricingStabilizationDynamic() {}

void PricingStabilizationDynamic::solve(double  _eta_out,
                                        double* _pi_out,
                                        double* _lhs) {
    k = 0.0;
    double result_sep;
    bool   mispricing = true;
    update = 0;

    std::copy(_pi_out,
              _pi_out + solver->reformulation_model.get_nb_constraints(),
              pi_out.begin());
    eta_out = _eta_out;
    iterations++;
    update_stab_center = false;

    do {
        k += 1.0;
        alphabar = hasstabcenter ? CC_MAX(0.0, 1.0 - k * (1 - alpha)) : 0.0;
        compute_pi_eta_sep(alphabar);
        OptimalSolution<double> aux_sol =
            std::move(solver->pricing_algorithm(pi_sep.data()));
        result_sep = solver->compute_lagrange(aux_sol, pi_sep.data());
        reduced_cost =
            solver->compute_reduced_cost(aux_sol, pi_out.data(), _lhs);

        if (reduced_cost <= -1e-6) {
            solver->compute_subgradient(aux_sol, subgradient.data());
            adjust_alpha();
            sol = std::move(aux_sol);
            alpha = alphabar;
            update = 1;
            mispricing = false;
        }
    } while (mispricing && alphabar > 0.0);

    if (result_sep > eta_in) {
        hasstabcenter = 1;
        eta_in = result_sep;
        pi_in = pi_sep;
        update_stab_center = true;
        if (iterations % 10 == 0 && eta_out - 1e-2 <= solver->UB) {
            solver->calculate_constLB(pi_sep.data());
            solver->evaluate_nodes(pi_sep.data());
        }
    }

    if (iterations % solver->convex_constr_id == 0) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }
}

/**
@brief Construct a new Pricing Stabilization Hybrid:: Pricing Stabilization
Hybrid object
 *
 * @param pricer_solver
 */
PricingStabilizationHybrid::PricingStabilizationHybrid(
    PricerSolverBase* pricer_solver)
    : PricingStabilizationDynamic(pricer_solver),
      subgradient_in(pricer_solver->convex_constr_id + 1) {}

PricingStabilizationHybrid::~PricingStabilizationHybrid() {}

void PricingStabilizationHybrid::solve(double  _eta_out,
                                       double* _pi_out,
                                       double* _lhs) {
    update = 0;
    bool stabilized = false;
    std::copy(_pi_out,
              _pi_out + solver->reformulation_model.get_nb_constraints(),
              pi_out.begin());
    eta_out = _eta_out;
    iterations++;

    do {
        update_hybrid();

        stabilized = is_stabilized();

        for (int i = 0; i < solver->reformulation_model.get_nb_constraints();
             ++i) {
            pi_sep[i] = compute_dual(i);
        }

        OptimalSolution<double> aux_sol;
        aux_sol = solver->pricing_algorithm(pi_sep.data());

        eta_sep = solver->compute_lagrange(aux_sol, pi_sep.data());
        reduced_cost =
            solver->compute_reduced_cost(aux_sol, pi_out.data(), _lhs);

        update_stabcenter(aux_sol);

        if (reduced_cost < -1e-6) {
            if (in_mispricing_schedule) {
                in_mispricing_schedule = 0;
            }
            sol = std::move(aux_sol);
            update_subgradientproduct();
            update_alpha();
        } else {
            if (stabilized) {
                in_mispricing_schedule = 1;
                update_alpha_misprice();
            } else {
                in_mispricing_schedule = 0;
            }
        }
    } while (in_mispricing_schedule && stabilized);

    if (iterations % solver->convex_constr_id == 0) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }
}

extern "C" {
PricingStabilizationBase* new_pricing_stabilization(PricerSolver* solver,
                                                    Parms*        parms) {
    switch (parms->stab_technique) {
        case stab_wentgnes:
            return new PricingStabilizationStat(solver);
        case stab_dynamic:
            return new PricingStabilizationDynamic(solver);
        case stab_hybrid:
            return new PricingStabilizationHybrid(solver);
        case no_stab:
            return new PricingStabilizationBase(solver);
        default:
            return new PricingStabilizationStat(solver);
    }
}

void delete_pricing_stabilization(PricingStabilization* pricing_stab_solver) {
    if (pricing_stab_solver) {
        delete pricing_stab_solver;
    }
}

double call_get_reduced_cost(PricingStabilizationBase* p) {
    return p->get_reduced_cost();
}

double call_get_eta_in(PricingStabilizationBase* solver) {
    return solver->get_eta_in();
}

int call_stopping_criteria(PricingStabilizationBase* solver) {
    return solver->stopping_criteria();
}

int call_get_update_stab_center(PricingStabilizationBase* solver) {
    return solver->get_update_stab_center();
}
}