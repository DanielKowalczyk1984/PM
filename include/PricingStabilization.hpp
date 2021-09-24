#ifndef __PRICINGSTABILIZATION_H__
#define __PRICINGSTABILIZATION_H__

#include <cstddef>              // for size_t
#include <memory>               // for unique_ptr
#include <vector>               // for vector
#include "PricingSolution.hpp"  // for PricingSolution
struct PricerSolverBase;        // lines 9-9
class PricingStabilizationBase {
   public:
    explicit PricingStabilizationBase(PricerSolverBase*        _solver,
                                      std::span<const double>& _pi_out);
    PricingStabilizationBase(PricingStabilizationBase&&) = default;
    PricingStabilizationBase(const PricingStabilizationBase&) = delete;
    PricingStabilizationBase& operator=(PricingStabilizationBase&&) = delete;
    PricingStabilizationBase& operator=(const PricingStabilizationBase&) =
        delete;
    virtual ~PricingStabilizationBase() = default;
    virtual std::unique_ptr<PricingStabilizationBase> clone(
        PricerSolverBase*,
        std::span<const double>&);

    static constexpr double EPS_RC = -1e-10;
    static constexpr double ETA_DIFF = 1e-4;
    static constexpr double EPS_STAB = 1e-9;
    static constexpr int    RC_FIXING_RATE = 50;

   public:
    PricerSolverBase* solver;
    PricingSolution<> sol;

    double reduced_cost{};
    double eta_in{};
    double eta_sep{};
    double eta_out{};

    // std::vector<double>& pi_out;
    std::span<const double>& pi_out;
    std::vector<double>      pi_in;
    std::vector<double>      pi_sep;

    int    update_stab_center{};
    size_t iterations{};
    size_t previous_fix{};

    bool continueLP{};

    PricingSolution<>& get_sol();
    bool               get_update_stab_center();
    double             get_reduced_cost();
    double             get_eps_stab_solver();
    virtual double     get_eta_in();
    virtual double     get_eta_sep();
    virtual void       solve(double eta_out, double* _lhs);
    virtual bool       stopping_criteria();
    virtual void       update_duals();
    virtual bool       reduced_cost_fixing();
    virtual void       remove_constraints(int first, int nb_del);
    virtual void       update_continueLP(int _continueLP);
    virtual int        do_reduced_cost_fixing();
    virtual void       update_continueLP(double obj);
    virtual void       set_alpha([[maybe_unused]] double _alpha){};
};

class PricingStabilizationStat : public PricingStabilizationBase {
   public:
    explicit PricingStabilizationStat(PricerSolverBase*        _solver,
                                      std::span<const double>& _pi_out);
    PricingStabilizationStat(PricingStabilizationStat&&) = default;
    PricingStabilizationStat(const PricingStabilizationStat&) = delete;
    PricingStabilizationStat& operator=(PricingStabilizationStat&&) = delete;
    PricingStabilizationStat& operator=(const PricingStabilizationStat&) =
        delete;
    ~PricingStabilizationStat() override = default;

    static constexpr double ALPHA = 0.8;
    static constexpr double ALPHA_CHG = 0.1;
    static constexpr double ALPHA_MAX = 0.99;

    double alpha{ALPHA};
    double alphabar{};
    int    update{};
    int    k{};
    int    hasstabcenter{};
    void   compute_pi_eta_sep(double _alpha_bar);

   public:
    void solve(double _eta_out, double* _lhs_coeff) override;

    double get_eta_in() final;
    bool   stopping_criteria() final;
    void   update_duals() override;
    void   remove_constraints(int first, int nb_del) override;
    void   set_alpha(double _alpha) override;
    std::unique_ptr<PricingStabilizationBase> clone(
        PricerSolverBase*,
        std::span<const double>&) override;
};

class PricingStabilizationDynamic : public PricingStabilizationStat {
   public:
    explicit PricingStabilizationDynamic(PricerSolverBase*        _solver,
                                         std::span<const double>& _pi_out);
    PricingStabilizationDynamic(PricingStabilizationDynamic&&) = default;
    PricingStabilizationDynamic(const PricingStabilizationDynamic&) = delete;
    PricingStabilizationDynamic& operator=(PricingStabilizationDynamic&&) =
        delete;
    PricingStabilizationDynamic& operator=(const PricingStabilizationDynamic&) =
        delete;
    ~PricingStabilizationDynamic() override = default;

    std::unique_ptr<PricingStabilizationBase> clone(
        PricerSolverBase*,
        std::span<const double>&) override;

    std::vector<double> subgradient;
    double              subgradientnorm{};
    double              subgradientproduct{};
    double              dualdiffnorm{};

    void solve(double _eta_out, double* _lhs_coeff) override;
    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

    void compute_subgradient(const PricingSolution<double>& _sol);
    void adjust_alpha();
};

class PricingStabilizationHybrid : public PricingStabilizationDynamic {
   public:
    explicit PricingStabilizationHybrid(PricerSolverBase*        pricer_solver,
                                        std::span<const double>& _pi_out);

    PricingStabilizationHybrid(PricingStabilizationHybrid&&) = default;
    PricingStabilizationHybrid(const PricingStabilizationHybrid&) = delete;
    PricingStabilizationHybrid& operator=(PricingStabilizationHybrid&&) =
        delete;
    PricingStabilizationHybrid& operator=(const PricingStabilizationHybrid&) =
        delete;
    ~PricingStabilizationHybrid() override = default;

    std::unique_ptr<PricingStabilizationBase> clone(
        PricerSolverBase*,
        std::span<const double>&) override;

    std::vector<double> subgradient_in;

    double beta{};
    double hybridfactor{};
    int    in_mispricing_schedule{};

    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

   private:
    void solve(double _eta_out, double* _lhs) override;

    void update_hybrid() {
        if (hasstabcenter && !in_mispricing_schedule && alpha > 0.0) {
            calculate_dualdiffnorm();
            calculate_beta();
            calculate_hybridfactor();
        }
    }

    void calculate_dualdiffnorm();

    void calculate_beta();

    void calculate_hybridfactor();

    bool is_stabilized();

    double compute_dual(auto i);

    void update_stabcenter(const PricingSolution<double>& _sol);

    void compute_subgradient_norm(const PricingSolution<double>& _sol);

    void update_subgradientproduct();

    void update_alpha_misprice();

    void update_alpha();
};

#endif  // __PRICINGSTABILIZATION_H__