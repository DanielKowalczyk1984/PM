#ifndef __PRICINGSTABILIZATION_H__
#define __PRICINGSTABILIZATION_H__

#include <algorithm>
#include <vector>
#include "OptimalSolution.hpp"
#include "PricerSolverBase.hpp"

class PricingStabilizationBase {
   public:
    PricingStabilizationBase(PricerSolverBase* _solver);
    PricingStabilizationBase(PricingStabilizationBase&&) = default;
    PricingStabilizationBase(const PricingStabilizationBase&) = default;
    PricingStabilizationBase& operator=(PricingStabilizationBase&&) = default;
    PricingStabilizationBase& operator=(const PricingStabilizationBase&) =
        default;
    ~PricingStabilizationBase();

   public:
    PricerSolverBase* solver;
    OptimalSolution<> sol;
    double            reduced_cost{};

    OptimalSolution<>& get_sol();
    double             get_reduced_cost();
    int                update_stab_center{};

    bool           get_update_stab_center();
    virtual double get_eta_in();
    virtual void   solve(double eta_out, double* _pi_out, double* _lhs);
    virtual int    stopping_criteria();
};

class PricingStabilizationStat : public PricingStabilizationBase {
   public:
    PricingStabilizationStat(PricerSolverBase* _solver);
    PricingStabilizationStat(PricingStabilizationStat&&) = default;
    PricingStabilizationStat(const PricingStabilizationStat&) = default;
    PricingStabilizationStat& operator=(PricingStabilizationStat&&) = default;
    PricingStabilizationStat& operator=(const PricingStabilizationStat&) =
        default;
    ~PricingStabilizationStat();

    std::vector<double> pi_in;
    std::vector<double> pi_out;
    std::vector<double> pi_sep;

    double alpha{0.8};
    double alphabar{};
    double eta_in{};
    double eta_out{};
    double eta_sep{};

    int iterations{};
    int update{};

    int k{};
    // int node_stab{};
    int hasstabcenter{};
    // int in_mispricing_schedule{};
    // int update_stab_center{};

    // PricerSolverBase* solver;
    // OptimalSolution<> sol{};

    int nb_rows;

    void compute_pi_eta_sep(double _alpha_bar) {
        int    i;
        double beta_bar = 1.0 - _alpha_bar;

        for (i = 0; i < solver->nb_jobs + 1; ++i) {
            pi_sep[i] = _alpha_bar * pi_in[i] + (1.0 - _alpha_bar) * pi_out[i];
        }

        eta_sep = _alpha_bar * (eta_in) + beta_bar * (eta_out);
    }

   public:
    virtual void solve(double  _eta_out,
                       double* _pi_out,
                       double* _lhs_coeff) override;

    virtual double get_eta_in() final;
    virtual int    stopping_criteria() final;
    // double diff_eta_in_eta_out();
};

class PricingStabilizationDynamic : public PricingStabilizationStat {
   public:
    PricingStabilizationDynamic(PricerSolverBase* _solver);
    PricingStabilizationDynamic(PricingStabilizationDynamic&&) = default;
    PricingStabilizationDynamic(const PricingStabilizationDynamic&) = default;
    PricingStabilizationDynamic& operator=(PricingStabilizationDynamic&&) =
        default;
    PricingStabilizationDynamic& operator=(const PricingStabilizationDynamic&) =
        default;
    ~PricingStabilizationDynamic();
    std::vector<double> subgradient;
    double              subgradientnorm{};
    double              subgradientproduct{};
    double              dualdiffnorm{};

    virtual void solve(double  _eta_out,
                       double* _pi_out,
                       double* _lhs_coeff) override;

    void compute_subgradient(const OptimalSolution<double>& sol) {
        solver->compute_subgradient(sol, subgradient.data());

        // subgradientnorm = 0.0;

        // for (int i = 0; i < nb_rows; ++i) {
        //     double sqr = SQR(subgradient_in[i]);

        //     if (sqr > 0.00001) {
        //         subgradientnorm += sqr;
        //     }
        // }

        // subgradientnorm = SQRT(subgradientnorm);
    }

    void adjust_alpha() {
        double sum = 0.0;

        for (int i = 0; i < solver->reformulation_model.get_nb_constraints();
             ++i) {
            sum += subgradient[i] * (pi_out[i] - pi_in[i]);
        }

        if (sum > 0) {
            alphabar = std::max<double>(0, alphabar - 0.1);
        } else {
            alphabar = std::min<double>(0.9, alphabar + (1 - alphabar) * 0.05);
        }
    }

    // double hybridfactor{};
    // double beta{};

    // void update_alpha() {
    //     if (subgradientproduct > 0.0) {
    //         alpha = std::max<double>(0, alpha - 0.1);
    //     } else {
    //         alpha = std::min<double>(0.99, alpha + (1 - alpha) * 0.1);
    //     }
    // }

    // void update_alpha_misprice() {
    //     k++;
    //     alphabar = std::max<double>(0.0, 1 - k * (1.0 - alpha));
    // }

    // int is_stabilized() {
    //     if (in_mispricing_schedule) {
    //         return alphabar > 0.0;
    //     }
    //     return alpha > 0.0;
    // }

    // int calculate_dualdiffnorm() {
    //     int val = 0;

    //     dualdiffnorm = 0.0;

    //     for (int i = 0; i < nb_rows; ++i) {
    //         double dualdiff = SQR(pi_in[i] - pi_out[i]);
    //         if (dualdiff > 0.00001) {
    //             dualdiffnorm += dualdiff;
    //         }
    //     }

    //     dualdiffnorm = SQRT(dualdiffnorm);

    //     return val;
    // }

    // int calculate_beta() {
    //     int val = 0;

    //     beta = 0.0;
    //     for (int i = 0; i < nb_rows; ++i) {
    //         double dualdiff = ABS(pi_out[i] - pi_in[i]);
    //         double product = dualdiff * ABS(subgradient_in[i]);

    //         if (product > 0.000001) {
    //             beta += product;
    //         }
    //     }

    //     if (subgradientnorm > 0.00001) {
    //         beta = beta / (subgradientnorm * dualdiffnorm);
    //     }

    //     return val;
    // }

    // int calculate_hybridfactor() {
    //     int val = 0;

    //     double aux_norm = 0.0;
    //     for (int i = 0; i < nb_rows; ++i) {
    //         double aux_double = SQR(
    //             (beta - 1.0) * (pi_out[i] - pi_in[i]) +
    //             beta * (subgradient_in[i] * dualdiffnorm / subgradientnorm));
    //         if (aux_double > 0.00001) {
    //             aux_norm += aux_double;
    //         }
    //     }
    //     aux_norm = SQRT(aux_norm);

    //     hybridfactor = ((1 - alpha) * dualdiffnorm) / aux_norm;

    //     return val;
    // }

    // double compute_dual(int i) {
    //     double usedalpha = alpha;
    //     double usedbeta = beta;

    //     if (in_mispricing_schedule) {
    //         usedalpha = alphabar;
    //         usedbeta = 0.0;
    //     }

    //     if (hasstabcenter && (usedbeta == 0.0 || usedalpha == 0.0)) {
    //         return usedalpha * pi_in[i] + (1.0 - usedalpha) * pi_out[i];
    //     } else if (hasstabcenter && usedbeta > 0.0) {
    //         double dual =
    //             pi_in[i] +
    //             hybridfactor *
    //                 (beta * (pi_in[i] + subgradient_in[i] * dualdiffnorm /
    //                                         subgradientnorm) +
    //                  (1.0 - beta) * pi_out[i] - pi_in[i]);
    //         return CC_MAX(dual, 0.0);
    //     }

    //     return pi_out[i];
    // }

    // int row_getDual(int i) {
    //     int val = 0;
    //     assert(i < nb_rows);

    //     pi_sep[i] = compute_dual(i);

    //     return val;
    // }

    // int update_subgradientproduct() {
    //     int val = 0;

    //     subgradientproduct = 0.0;
    //     for (int i = 0; i < nb_rows; ++i) {
    //         subgradientproduct -= (pi_out[i] - pi_in[i]) * subgradient_in[i];
    //     }

    //     return val;
    // }

    // int update_stabcenter(const OptimalSolution<double>& sol) {
    //     int val = 0;

    //     if (eta_sep > eta_in) {
    //         pi_in = pi_sep;
    //         compute_subgradient(sol);
    //         eta_in = eta_sep;
    //         hasstabcenter = 1;
    //     }

    //     return val;
    // }
};

class PricingStabilizationHybrid : public PricingStabilizationDynamic {
   public:
    PricingStabilizationHybrid(PricerSolverBase* pricer_solver);
    PricingStabilizationHybrid(PricingStabilizationHybrid&&) = default;
    PricingStabilizationHybrid(const PricingStabilizationHybrid&) = default;
    PricingStabilizationHybrid& operator=(PricingStabilizationHybrid&&) =
        default;
    PricingStabilizationHybrid& operator=(const PricingStabilizationHybrid&) =
        default;
    ~PricingStabilizationHybrid();

    std::vector<double> subgradient_in;

    double beta{};
    double hybridfactor{};
    int    in_mispricing_schedule{};

   private:
    virtual void solve(double _eta_out, double* _pi_out, double* _lhs) override;

    void update_hybrid() {
        if (hasstabcenter && !in_mispricing_schedule && alpha > 0.0) {
            calculate_dualdiffnorm();
            calculate_beta();
            calculate_hybridfactor();
        }
    }

    void calculate_dualdiffnorm() {
        dualdiffnorm = 0.0;

        for (int i = 0; i < solver->reformulation_model.get_nb_constraints();
             ++i) {
            double dualdiff = SQR(pi_in[i] - pi_out[i]);
            if (dualdiff > 0.00001) {
                dualdiffnorm += dualdiff;
            }
        }

        dualdiffnorm = SQRT(dualdiffnorm);
    }

    void calculate_beta() {
        beta = 0.0;
        for (int i = 0; i < solver->nb_jobs; ++i) {
            double dualdiff = ABS(pi_out[i] - pi_in[i]);
            double product = dualdiff * std::abs(subgradient_in[i]);

            if (product > 1e-6) {
                beta += product;
            }
        }

        if (subgradientnorm > 1e-5) {
            beta = beta / (subgradientnorm * dualdiffnorm);
        }
    }

    void calculate_hybridfactor() {
        double aux_norm = 0.0;
        for (int i = 0; i < solver->nb_jobs; ++i) {
            double aux_double = SQR(
                (beta - 1.0) * (pi_out[i] - pi_in[i]) +
                beta * (subgradient_in[i] * dualdiffnorm / subgradientnorm));
            if (aux_double > 1e-6) {
                aux_norm += aux_double;
            }
        }
        aux_norm = SQRT(aux_norm);

        hybridfactor = ((1 - alpha) * dualdiffnorm) / aux_norm;
    }

    bool is_stabilized() {
        if (in_mispricing_schedule) {
            return alphabar > 0.0;
        }
        return alpha > 0.0;
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

    void update_stabcenter(const OptimalSolution<double>& sol) {
        if (eta_sep > eta_in) {
            pi_in = pi_sep;
            compute_subgradient_norm(sol);
            eta_in = eta_sep;
            hasstabcenter = 1;
            update_stab_center = true;
        }
    }

    void compute_subgradient_norm(const OptimalSolution<double>& sol) {
        // double* subgradient_in = &g_array_index(pd->subgradient_in, double,
        // 0); double* rhs = &g_array_index(pd->rhs, double, 0);
        // fill_dbl(subgradient_in, pd->nb_rows, 1.0);
        // subgradient_in[pd->nb_jobs] = 0.0;
        solver->compute_subgradient(sol, subgradient_in.data());

        // for (guint i = 0; i < sol.jobs->len; i++) {
        //     Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol.jobs,
        //     i)); subgradient_in[tmp_j->job] += rhs[pd->nb_jobs] * 1.0;
        // }

        subgradientnorm = 0.0;

        for (int i = 0; i < solver->reformulation_model.get_nb_constraints();
             ++i) {
            double sqr = SQR(subgradient_in[i]);

            if (sqr > 0.00001) {
                subgradientnorm += sqr;
            }
        }

        subgradientnorm = SQRT(subgradientnorm);
    }
    void update_subgradientproduct() {
        subgradientproduct = 0.0;
        for (int i = 0; i < solver->reformulation_model.get_nb_constraints();
             ++i) {
            subgradientproduct += (pi_out[i] - pi_in[i]) * subgradient_in[i];
        }
    }

    void update_alpha_misprice() {
        k++;
        alphabar = CC_MAX(0.0, 1 - k * (1.0 - alpha));
    }

    void update_alpha() {
        if (subgradientproduct > 0.0) {
            alpha = std::max(0.0, alpha - 0.1);
        } else {
            alpha = std::min(0.99, alpha + (1.0 - alpha) * 0.1);
        }
    }
};

#endif  // __PRICINGSTABILIZATION_H__