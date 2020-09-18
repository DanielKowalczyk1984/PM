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
    : solver(_solver),
      pi_in(_solver->convex_constr_id + 1, 0.0),
      pi_sep(_solver->convex_constr_id + 1, 0.0) {}

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
    std::copy(_pi_out,
              _pi_out + solver->reformulation_model.get_nb_constraints(),
              pi_sep.begin());
    solver->calculate_constLB(pi_sep.data());
    sol = std::move(solver->pricing_algorithm(pi_sep.data()));
    reduced_cost = solver->compute_reduced_cost(sol, pi_sep.data(), _lhs);
    eta_sep = solver->compute_lagrange(sol, pi_sep.data());
    update_stab_center = false;
    eta_out = _eta_out;
    continueLP = (eta_in < eta_out);
    if (eta_sep > eta_in) {
        eta_in = eta_sep;
        pi_in = pi_sep;
        update_stab_center = true;
    }
    iterations++;
    // if (iterations % 20 == 0 && _eta_out - 1e-2 <= solver->UB) {
    // solver->calculate_constLB(_pi_out);
    // solver->evaluate_nodes(_pi_out);
    // }
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
    return eta_in;
}

double PricingStabilizationBase::get_eta_sep() {
    return eta_sep;
}

int PricingStabilizationBase::stopping_criteria() {
    return reduced_cost < -1e-6;
}

void PricingStabilizationBase::reduced_cost_fixing() {
    previous_fix = iterations;
    solver->calculate_constLB(pi_sep.data());
    solver->evaluate_nodes(pi_sep.data());
}

void PricingStabilizationBase::remove_constraints(int first, int nb_del) {
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
}

void PricingStabilizationBase::update_continueLP(int _continueLP) {
    if (_continueLP) {
        continueLP = true;
    } else {
        continueLP = false;
    }
}

int PricingStabilizationBase::do_reduced_cost_fixing() {
    if (iterations == 1 || !continueLP || previous_fix <= iterations - 20) {
        fmt::print("reduced cost fixing test\n");
        return 1;
    }
    return 0;
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
      pi_out(_solver->convex_constr_id + 1, 0.0) {}

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
    bool mispricing = true;
    update = 0;

    std::copy(_pi_out,
              _pi_out + solver->reformulation_model.get_nb_constraints(),
              pi_out.begin());
    eta_out = _eta_out;
    iterations++;
    update_stab_center = false;

    do {
        k += 1.0;
        alphabar = (hasstabcenter && !(eta_in < eta_out - 0.02))
                       ? CC_MAX(0, 1.0 - k * (1.0 - alpha))
                       : 0.0;
        compute_pi_eta_sep(alphabar);
        solver->calculate_constLB(pi_sep.data());
        OptimalSolution<> aux_sol =
            std::move(solver->pricing_algorithm(pi_sep.data()));

        eta_sep = solver->compute_lagrange(aux_sol, pi_sep.data());
        reduced_cost =
            solver->compute_reduced_cost(aux_sol, pi_out.data(), _lhs_coeff);

        continueLP = (eta_sep < eta_out - 1e-2);

        if (reduced_cost < -1e-6) {
            sol = std::move(aux_sol);
            update = 1;
            mispricing = false;
        }
    } while (mispricing && alphabar > 0); /** mispricing check */

    if (eta_sep > eta_in) {
        hasstabcenter = 1;
        eta_in = eta_sep;
        pi_in = pi_sep;
        update_stab_center = true;
        // if (iterations % 50 == 0 && eta_out - 1e-2 <= solver->UB) {
        //     solver->calculate_constLB(pi_sep.data());
        //     solver->evaluate_nodes(pi_sep.data());
        // }
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

void PricingStabilizationStat::update_duals() {
    if (solver->reformulation_model.get_nb_constraints() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        pi_out.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        pi_in.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
    }
}

void PricingStabilizationStat::remove_constraints(int first, int nb_del) {
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
    auto it_out = pi_out.begin() + first;
    pi_out.erase(it_out, it_out + nb_del);
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
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
        auto aux_sol = std::move(solver->pricing_algorithm(pi_sep.data()));
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

void PricingStabilizationDynamic::update_duals() {
    if (solver->reformulation_model.get_nb_constraints() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        pi_out.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        pi_in.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        subgradient.resize(solver->reformulation_model.get_nb_constraints(),
                           0.0);
    }
}

void PricingStabilizationDynamic::remove_constraints(int first, int nb_del) {
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
    auto it_out = pi_out.begin() + first;
    pi_out.erase(it_out, it_out + nb_del);
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
    auto it_sub = subgradient.begin() + first;
    subgradient.erase(it_sub, it_sub + nb_del);
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

void PricingStabilizationHybrid::update_duals() {
    if (solver->reformulation_model.get_nb_constraints() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        pi_out.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        pi_in.resize(solver->reformulation_model.get_nb_constraints(), 0.0);
        subgradient_in.resize(solver->reformulation_model.get_nb_constraints(),
                              0.0);
        subgradient.resize(solver->reformulation_model.get_nb_constraints(),
                           0.0);
    }
}

void PricingStabilizationHybrid::remove_constraints(int first, int nb_del) {
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
    auto it_out = pi_out.begin() + first;
    pi_out.erase(it_out, it_out + nb_del);
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
    auto it_sub = subgradient.begin() + first;
    subgradient.erase(it_sub, it_sub + nb_del);
    auto it_sub_in = subgradient_in.begin() + first;
    subgradient_in.erase(it_sub_in, it_sub_in + nb_del);
}

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

void call_update_duals(PricingStabilizationBase* solver) {
    solver->update_duals();
}

void call_remove_constraints(PricingStabilizationBase* solver,
                             int                       first,
                             int                       nb_del) {
    solver->remove_constraints(first, nb_del);
}

int call_do_reduced_fixing(PricingStabilizationBase* solver) {
    return solver->do_reduced_cost_fixing();
}

void call_reduced_cost_fixing(PricingStabilizationBase* solver) {
    solver->reduced_cost_fixing();
}

void call_update_continueLP(PricingStabilizationBase* solver, double _eta_out) {
    solver->eta_out = _eta_out;
    solver->continueLP = (solver->eta_in < solver->eta_out);
}

int call_get_continueLP(PricingStabilizationBase* solver) {
    return solver->continueLP;
}

double call_get_eta_sep(PricingStabilizationBase* solver) {
    return solver->get_eta_sep();
}
}