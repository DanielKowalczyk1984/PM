#ifndef __PRICINGSTABILIZATION_H__
#define __PRICINGSTABILIZATION_H__

#include <algorithm>
#include <memory>
#include <vector>
#include "OptimalSolution.hpp"
// #include "OptimalSolution.hpp"
// #include "PricerSolverBase.hpp"
class PricerSolverBase;
class PricingStabilizationBase {
   public:
    explicit PricingStabilizationBase(PricerSolverBase* _solver);
    PricingStabilizationBase(PricingStabilizationBase&&) = default;
    PricingStabilizationBase(const PricingStabilizationBase&) = delete;
    PricingStabilizationBase& operator=(PricingStabilizationBase&&) = default;
    PricingStabilizationBase& operator=(const PricingStabilizationBase&) =
        delete;
    virtual ~PricingStabilizationBase();
    virtual std::unique_ptr<PricingStabilizationBase> clone(PricerSolverBase*);

    static constexpr double EPS_RC = -1e-10;
    static constexpr double ETA_DIFF = 1e-4;
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
    virtual void       solve(double                     eta_out,
                             const std::vector<double>& _pi_out,
                             double*                    _lhs);
    virtual int        stopping_criteria();
    virtual void       update_duals();
    virtual void       reduced_cost_fixing();
    virtual void       remove_constraints(int first, int nb_del);
    virtual void       update_continueLP(int _continueLP);
    virtual int        do_reduced_cost_fixing();
    virtual void       update_continueLP(double obj);
};

class PricingStabilizationStat : public PricingStabilizationBase {
   public:
    explicit PricingStabilizationStat(PricerSolverBase* _solver);
    PricingStabilizationStat(PricingStabilizationStat&&) = default;
    PricingStabilizationStat(const PricingStabilizationStat&) = delete;
    PricingStabilizationStat& operator=(PricingStabilizationStat&&) = default;
    PricingStabilizationStat& operator=(const PricingStabilizationStat&) =
        delete;
    ~PricingStabilizationStat() override;

    std::vector<double> pi_out;

    static constexpr double ALPHA = 0.8;
    static constexpr double ALPHA_CHG = 0.1;
    static constexpr double ALPHA_MAX = 0.99;

    double alpha{ALPHA};
    double alphabar{};

    int update{};

    int k{};
    int hasstabcenter{};

    int nb_rows;

    void compute_pi_eta_sep(double _alpha_bar);

   public:
    void solve(double                     _eta_out,
               const std::vector<double>& _pi_out,
               double*                    _lhs_coeff) override;

    double get_eta_in() final;
    int    stopping_criteria() final;
    void   update_duals() override;
    void   remove_constraints(int first, int nb_del) override;
    std::unique_ptr<PricingStabilizationBase> clone(PricerSolverBase*) override;
    // double diff_eta_in_eta_out();
};

class PricingStabilizationDynamic : public PricingStabilizationStat {
   public:
    explicit PricingStabilizationDynamic(PricerSolverBase* _solver);
    PricingStabilizationDynamic(PricingStabilizationDynamic&&) = default;
    PricingStabilizationDynamic(const PricingStabilizationDynamic&) = delete;
    PricingStabilizationDynamic& operator=(PricingStabilizationDynamic&&) =
        default;
    PricingStabilizationDynamic& operator=(const PricingStabilizationDynamic&) =
        delete;
    ~PricingStabilizationDynamic() override;
    std::unique_ptr<PricingStabilizationBase> clone(PricerSolverBase*) override;

    std::vector<double> subgradient;
    double              subgradientnorm{};
    double              subgradientproduct{};
    double              dualdiffnorm{};

    void solve(double                     _eta_out,
               const std::vector<double>& _pi_out,
               double*                    _lhs_coeff) override;
    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

    void compute_subgradient(const OptimalSolution<double>& _sol);
    void adjust_alpha();
};

class PricingStabilizationHybrid : public PricingStabilizationDynamic {
   public:
    explicit PricingStabilizationHybrid(PricerSolverBase* pricer_solver);
    PricingStabilizationHybrid(PricingStabilizationHybrid&&) = default;
    PricingStabilizationHybrid(const PricingStabilizationHybrid&) = delete;
    PricingStabilizationHybrid& operator=(PricingStabilizationHybrid&&) =
        default;
    PricingStabilizationHybrid& operator=(const PricingStabilizationHybrid&) =
        delete;
    ~PricingStabilizationHybrid() override;
    std::unique_ptr<PricingStabilizationBase> clone(PricerSolverBase*) override;

    std::vector<double> subgradient_in;

    double beta{};
    double hybridfactor{};
    int    in_mispricing_schedule{};

    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

   private:
    void solve(double                     _eta_out,
               const std::vector<double>& _pi_out,
               double*                    _lhs) override;

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

    void update_stabcenter(const OptimalSolution<double>& _sol);

    void compute_subgradient_norm(const OptimalSolution<double>& _sol);

    void update_subgradientproduct();

    void update_alpha_misprice();

    void update_alpha();
};

#endif  // __PRICINGSTABILIZATION_H__