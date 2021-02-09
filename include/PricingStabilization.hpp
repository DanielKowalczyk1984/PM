#ifndef __PRICINGSTABILIZATION_H__
#define __PRICINGSTABILIZATION_H__

#include <algorithm>
#include <vector>
#include "OptimalSolution.hpp"
#include "PricerSolverBase.hpp"

class PricingStabilizationBase {
   public:
    explicit PricingStabilizationBase(PricerSolverBase* _solver);
    PricingStabilizationBase(PricingStabilizationBase&&) = default;
    PricingStabilizationBase(const PricingStabilizationBase&) = default;
    PricingStabilizationBase& operator=(PricingStabilizationBase&&) = default;
    PricingStabilizationBase& operator=(const PricingStabilizationBase&) =
        default;
    virtual ~PricingStabilizationBase();

    static constexpr double EPS_RC = -1e-10;
    static constexpr double ETA_DIFF = 1e-6;
    static constexpr double EPS_STAB = 1e-9;
    static constexpr int    RC_FIXING_RATE = 20;

   public:
    PricerSolverBase* solver;
    OptimalSolution<> sol;

    double reduced_cost{};
    double eta_in{};
    double eta_sep{};
    double eta_out{};

    std::vector<double> pi_in;
    std::vector<double> pi_sep;

    int update_stab_center{};
    int iterations{};
    int previous_fix{};

    bool continueLP{};

    OptimalSolution<>& get_sol();
    bool               get_update_stab_center();
    double             get_reduced_cost();
    double             get_eps_stab_solver();
    virtual double     get_eta_in();
    virtual double     get_eta_sep();
    virtual void       solve(double eta_out, double* _pi_out, double* _lhs);
    virtual int        stopping_criteria();
    virtual void       update_duals();
    virtual void       reduced_cost_fixing();
    virtual void       remove_constraints(int first, int nb_del);
    virtual void       update_continueLP(int _continueLP);
    virtual int        do_reduced_cost_fixing();
};

class PricingStabilizationStat : public PricingStabilizationBase {
   public:
    explicit PricingStabilizationStat(PricerSolverBase* _solver);
    PricingStabilizationStat(PricingStabilizationStat&&) = default;
    PricingStabilizationStat(const PricingStabilizationStat&) = default;
    PricingStabilizationStat& operator=(PricingStabilizationStat&&) = default;
    PricingStabilizationStat& operator=(const PricingStabilizationStat&) =
        default;
    ~PricingStabilizationStat() override;

    std::vector<double> pi_out;

    static constexpr double ALPHA = 0.8;
    static constexpr double ALPHA_CHG = 0.1;
    static constexpr double ALPHA_MAX = 0.99;

    double alpha{ALPHA};
    double alphabar{};

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
        double beta_bar = 1.0 - _alpha_bar;

        for (int i = 0; i < solver->reformulation_model.get_nb_constraints();
             ++i) {
            pi_sep[i] = _alpha_bar * pi_in[i] + (1.0 - _alpha_bar) * pi_out[i];
        }

        eta_sep = _alpha_bar * (eta_in) + beta_bar * (eta_out);
    }

   public:
    void solve(double _eta_out, double* _pi_out, double* _lhs_coeff) override;

    double get_eta_in() final;
    int    stopping_criteria() final;
    void   update_duals() override;
    void   remove_constraints(int first, int nb_del) override;
    // double diff_eta_in_eta_out();
};

class PricingStabilizationDynamic : public PricingStabilizationStat {
   public:
    explicit PricingStabilizationDynamic(PricerSolverBase* _solver);
    PricingStabilizationDynamic(PricingStabilizationDynamic&&) = default;
    PricingStabilizationDynamic(const PricingStabilizationDynamic&) = default;
    PricingStabilizationDynamic& operator=(PricingStabilizationDynamic&&) =
        default;
    PricingStabilizationDynamic& operator=(const PricingStabilizationDynamic&) =
        default;
    ~PricingStabilizationDynamic() override;
    std::vector<double> subgradient;
    double              subgradientnorm{};
    double              subgradientproduct{};
    double              dualdiffnorm{};

    void solve(double _eta_out, double* _pi_out, double* _lhs_coeff) override;
    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

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
            alphabar = std::max<double>(0, alphabar - ALPHA_CHG);
        } else {
            alphabar = std::min<double>(ALPHA_MAX,
                                        alphabar + (1 - alphabar) * ALPHA_CHG);
        }
    }
};

class PricingStabilizationHybrid : public PricingStabilizationDynamic {
   public:
    explicit PricingStabilizationHybrid(PricerSolverBase* pricer_solver);
    PricingStabilizationHybrid(PricingStabilizationHybrid&&) = default;
    PricingStabilizationHybrid(const PricingStabilizationHybrid&) = default;
    PricingStabilizationHybrid& operator=(PricingStabilizationHybrid&&) =
        default;
    PricingStabilizationHybrid& operator=(const PricingStabilizationHybrid&) =
        default;
    ~PricingStabilizationHybrid() override;

    std::vector<double> subgradient_in;

    double beta{};
    double hybridfactor{};
    int    in_mispricing_schedule{};

    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

   private:
    void solve(double _eta_out, double* _pi_out, double* _lhs) override;

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
            if (dualdiff > EPS_STAB) {
                dualdiffnorm += dualdiff;
            }
        }

        dualdiffnorm = SQRT(dualdiffnorm);
    }

    void calculate_beta() {
        beta = 0.0;
        for (int i = 0; i < solver->convex_constr_id; ++i) {
            double dualdiff = ABS(pi_out[i] - pi_in[i]);
            double product = dualdiff * std::abs(subgradient_in[i]);

            if (product > EPS_STAB) {
                beta += product;
            }
        }

        if (subgradientnorm > EPS_STAB) {
            beta = beta / (subgradientnorm * dualdiffnorm);
        }
    }

    void calculate_hybridfactor() {
        double aux_norm = 0.0;
        for (int i = 0; i < solver->convex_constr_id; ++i) {
            double aux_double = SQR(
                (beta - 1.0) * (pi_out[i] - pi_in[i]) +
                beta * (subgradient_in[i] * dualdiffnorm / subgradientnorm));
            if (aux_double > EPS_STAB) {
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
            solver->calculate_constLB(pi_sep.data());
            solver->evaluate_nodes(pi_sep.data());
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

            if (sqr > EPS_STAB) {
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
            alpha = std::max(0.0, alpha - ALPHA_CHG);
        } else {
            alpha = std::min(ALPHA_MAX, alpha + (1.0 - alpha) * ALPHA_CHG);
        }
    }
};

#endif  // __PRICINGSTABILIZATION_H__