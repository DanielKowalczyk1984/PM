#ifndef __PRICINGSTABILIZATION_H__
#define __PRICINGSTABILIZATION_H__

#include <algorithm>
#include <vector>
#include "OptimalSolution.hpp"
#include "PricerSolverBase.hpp"
#include "fmt/core.h"

class PricingStabilization {
   public:
    PricingStabilization(PricerSolverBase* _solver);
    PricingStabilization(PricingStabilization&&) = default;
    PricingStabilization(const PricingStabilization&) = default;
    PricingStabilization& operator=(PricingStabilization&&) = default;
    PricingStabilization& operator=(const PricingStabilization&) = default;
    ~PricingStabilization();

   private:
    std::vector<double> pi_in;
    std::vector<double> subgradient;
    std::vector<double> pi_out;
    std::vector<double> pi_sep;
    std::vector<double> subgradient_in;

    double dualdiffnorm{};
    double hybridfactor{};
    double subgradientnorm{};
    double alpha{0.8};
    double alphabar{};
    double beta{};
    double eta_in{};
    double subgradientproduct{};
    double eta_out{};
    double eta_sep{};
    double reduced_cost{};

    int k{};
    int node_stab{};
    int hasstabcenter{};
    int in_mispricing_schedule{};
    int update_stab_center{};
    int update{};
    int iterations{};

    PricerSolverBase* solver;
    OptimalSolution<> sol{};

    int nb_rows;

    void update_alpha() {
        if (subgradientproduct > 0.0) {
            alpha = std::max<double>(0, alpha - 0.1);
        } else {
            alpha = std::min<double>(0.99, alpha + (1 - alpha) * 0.1);
        }
    }

    void update_alpha_misprice() {
        k++;
        alphabar = std::max<double>(0.0, 1 - k * (1.0 - alpha));
    }

    int is_stabilized() {
        if (in_mispricing_schedule) {
            return alphabar > 0.0;
        }
        return alpha > 0.0;
    }

    int calculate_dualdiffnorm() {
        int val = 0;

        dualdiffnorm = 0.0;

        for (int i = 0; i < nb_rows; ++i) {
            double dualdiff = SQR(pi_in[i] - pi_out[i]);
            if (dualdiff > 0.00001) {
                dualdiffnorm += dualdiff;
            }
        }

        dualdiffnorm = SQRT(dualdiffnorm);

        return val;
    }

    int calculate_beta() {
        int val = 0;

        beta = 0.0;
        for (int i = 0; i < nb_rows; ++i) {
            double dualdiff = ABS(pi_out[i] - pi_in[i]);
            double product = dualdiff * ABS(subgradient_in[i]);

            if (product > 0.000001) {
                beta += product;
            }
        }

        if (subgradientnorm > 0.00001) {
            beta = beta / (subgradientnorm * dualdiffnorm);
        }

        return val;
    }

    int calculate_hybridfactor() {
        int val = 0;

        double aux_norm = 0.0;
        for (int i = 0; i < nb_rows; ++i) {
            double aux_double = SQR(
                (beta - 1.0) * (pi_out[i] - pi_in[i]) +
                beta * (subgradient_in[i] * dualdiffnorm / subgradientnorm));
            if (aux_double > 0.00001) {
                aux_norm += aux_double;
            }
        }
        aux_norm = SQRT(aux_norm);

        hybridfactor = ((1 - alpha) * dualdiffnorm) / aux_norm;

        return val;
    }

    int update_hybrid() {
        int val = 0;

        if (hasstabcenter && !in_mispricing_schedule && alpha > 0.0) {
            calculate_dualdiffnorm();
            calculate_beta();
            calculate_hybridfactor();
        }

        return val;
    }

    int update_node() {
        int val = 0;
        // if (node_stab != id) {
        //     node_stab = id;
        //     k = 0;
        //     alpha = 0.0;
        //     hasstabcenter = 0;
        //     eta_in = 0.0;
        //     in_mispricing_schedule = 0;
        // }

        return val;
    }

    double compute_dual(int i) {
        double usedalpha = alpha;
        double usedbeta = beta;

        if (in_mispricing_schedule) {
            usedalpha = alphabar;
            usedbeta = 0.0;
        }

        if (hasstabcenter && (usedbeta == 0.0 || usedalpha == 0.0)) {
            return usedalpha * pi_in[i] + (1.0 - usedalpha) * pi_out[i];
        } else if (hasstabcenter && usedbeta > 0.0) {
            double dual =
                pi_in[i] +
                hybridfactor *
                    (beta * (pi_in[i] + subgradient_in[i] * dualdiffnorm /
                                            subgradientnorm) +
                     (1.0 - beta) * pi_out[i] - pi_in[i]);
            return CC_MAX(dual, 0.0);
        }

        return pi_out[i];
    }

    int row_getDual(int i) {
        int val = 0;
        assert(i < nb_rows);

        pi_sep[i] = compute_dual(i);

        return val;
    }

    void compute_subgradient(const OptimalSolution<double>& sol, double* rhs) {
        std::fill(subgradient_in.begin(), subgradient.end(), 1.0);
        subgradient_in[solver->nb_jobs] = 0.0;

        for (guint i = 0; i < sol.jobs->len; i++) {
            Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol.jobs, i));
            subgradient_in[tmp_j->job] += rhs[solver->nb_jobs] * 1.0;
        }

        subgradientnorm = 0.0;

        for (int i = 0; i < nb_rows; ++i) {
            double sqr = SQR(subgradient_in[i]);

            if (sqr > 0.00001) {
                subgradientnorm += sqr;
            }
        }

        subgradientnorm = SQRT(subgradientnorm);
    }

    int update_subgradientproduct() {
        int val = 0;

        subgradientproduct = 0.0;
        for (int i = 0; i < nb_rows; ++i) {
            subgradientproduct -= (pi_out[i] - pi_in[i]) * subgradient_in[i];
        }

        return val;
    }

    int update_stabcenter(const OptimalSolution<double>& sol) {
        int val = 0;

        if (eta_sep > eta_in) {
            pi_in = pi_sep;
            // compute_subgradient(sol);
            eta_in = eta_sep;
            hasstabcenter = 1;
        }

        return val;
    }

    void adjust_alpha() {
        double sum = 0.0;

        for (int i = 0; i < solver->nb_jobs; ++i) {
            sum += subgradient[i] * (pi_out[i] - pi_in[i]);
        }

        if (sum > 0) {
            alpha = std::max<double>(0, alpha - 0.1);
        } else {
            alpha = std::min<double>(0.9, alpha + (1 - alpha) * 0.05);
        }
    }

    void compute_pi_eta_sep(double _alpha_bar) {
        int    i;
        double beta_bar = 1.0 - _alpha_bar;

        for (i = 0; i < solver->nb_jobs + 1; ++i) {
            pi_sep[i] = _alpha_bar * pi_in[i] + (1.0 - _alpha_bar) * pi_out[i];
        }

        eta_sep = _alpha_bar * (eta_in) + beta_bar * (eta_out);
    }

   public:
    int solve_stab(double _eta_out, double* _pi_out, double* _lhs_coeff);

    OptimalSolution<>& get_sol();

    double get_reduced_cost();

    double get_eta_in();
    double get_diff_eta_in_eta_out();
};

#endif  // __PRICINGSTABILIZATION_H__