#include "PricingStabilization.hpp"
#include <cmath>

PricingStabilization::PricingStabilization(PricerSolverBase* _solver)
    : pi_in(_solver->nb_jobs + 1),
      pi_out(_solver->nb_jobs + 1),
      pi_sep(_solver->nb_jobs + 1),
      solver(_solver) {}

PricingStabilization::~PricingStabilization() {}

int PricingStabilization::solve_stab(double  _eta_out,
                                     double* _pi_out,
                                     double* _lhs_coeff) {
    int    val = 0;
    double k = 0.0;
    double alpha_bar;
    bool   mispricing = true;
    double result_sep;
    update = 0;

    std::copy(_pi_out, _pi_out + solver->nb_jobs + 1, pi_out.begin());
    eta_out = _eta_out;
    iterations++;

    do {
        k += 1.0;
        alpha_bar = hasstabcenter ? CC_MAX(0, 1.0 - k * (1.0 - alpha)) : 0.0;
        compute_pi_eta_sep(alpha_bar);
        OptimalSolution<> aux_sol =
            std::move(solver->pricing_algorithm(pi_sep.data()));

        result_sep = solver->compute_lagrange(aux_sol, pi_sep.data());
        reduced_cost =
            solver->compute_reduced_cost(aux_sol, pi_out.data(), _lhs_coeff);

        if (reduced_cost < -1e-6) {
            sol = std::move(aux_sol);
            CCcheck_val_2(val, "Failed in construct_sol_stab");
            update = 1;
            mispricing = false;
        }
    } while (mispricing && alpha_bar > 0); /** mispricing check */

    if (result_sep > eta_in) {
        hasstabcenter = 1;
        eta_in = result_sep;
        pi_in = pi_sep;
        update_stab_center = 1;
    } else {
        update_stab_center = 0;
    }

    if (iterations % solver->nb_jobs == 0) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }

CLEAN:
    return val;
}

OptimalSolution<>& PricingStabilization::get_sol() {
    return sol;
}

double PricingStabilization::get_reduced_cost() {
    return reduced_cost;
}

double PricingStabilization::get_eta_in() {
    return eta_in;
}

double PricingStabilization::get_diff_eta_in_eta_out() {
    return std::abs(eta_out - eta_in);
}

extern "C" double call_C_get_diff_eta(PricingStabilization* p) {
    return p->get_diff_eta_in_eta_out();
}

extern "C" PricingStabilization* new_pricing_stabilization(
    PricerSolver* solver) {
    return new PricingStabilization(solver);
}

extern "C" void delete_pricing_stabilization(
    PricingStabilization* pricing_stab_solver) {
    if (pricing_stab_solver) {
        delete pricing_stab_solver;
    }
}
