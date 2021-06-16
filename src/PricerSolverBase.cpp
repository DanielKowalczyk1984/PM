#include "PricerSolverBase.hpp"
#include <fmt/core.h>                                  // for print
#include <math.h>                                      // for fabs
#include <algorithm>                                   // for min, __fill_fn
#include <cstddef>                                     // for size_t
#include <ext/alloc_traits.h>                          // for __alloc_traits...
#include <limits>                                      // for numeric_limits
#include <memory>                                      // for __shared_ptr_a...
#include <range/v3/iterator/basic_iterator.hpp>        // for basic_iterator
#include <range/v3/iterator/unreachable_sentinel.hpp>  // for operator==
#include <range/v3/view/enumerate.hpp>                 // for enumerate_fn
#include <range/v3/view/view.hpp>                      // for operator|
#include <range/v3/view/zip.hpp>                       // for zip_view
#include <range/v3/view/zip_with.hpp>                  // for iter_zip_with_...
#include <span>                                        // for span
#include "Instance.h"                                  // for Instance
#include "Job.h"                                       // for Job
#include "gurobi_c++.h"                                // for GRBModel, GRBEnv
#include "gurobi_c.h"                                  // for GRB_INFEASIBLE

PricerSolverBase::PricerSolverBase(const Instance& instance)
    : jobs(instance.jobs),
      convex_constr_id(instance.nb_jobs),
      convex_rhs(instance.nb_machines),
      problem_name(),
      env(genv),
      model(*env),
      reformulation_model(instance.nb_jobs, instance.nb_machines),
      is_integer_solution(false),
      constLB(0.0),
      UB(std::numeric_limits<int>::max()),
      x_bar(std::vector<std::vector<double>>(
          instance.nb_jobs,
          std::vector<double>(instance.H_max + 1, 0.0))),
      z_bar(std::vector<std::vector<double>>(
          instance.nb_jobs,
          std::vector<double>(instance.H_max + 1, 0.0))) {
    try {
        model.set(GRB_IntParam_Method, GRB_METHOD_AUTO);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    } catch (const GRBException& e) {
        fmt::print("Error code = {}\n", e.getErrorCode());
        fmt::print(e.getMessage());
    } catch (...) {
        fmt::print("Exception during optimization\n");
    }
}

PricerSolverBase::PricerSolverBase(const PricerSolverBase& other)
    : jobs(other.jobs),
      convex_constr_id(other.convex_constr_id),
      convex_rhs(other.convex_rhs),
      problem_name(other.problem_name),
      env(other.env),
      model(*genv),
      reformulation_model(other.reformulation_model),
      is_integer_solution(other.is_integer_solution),
      constLB(other.constLB),
      UB(other.UB),
      x_bar(other.x_bar),
      z_bar(other.z_bar) {
    model.set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
}

PricerSolverBase::~PricerSolverBase() = default;

int PricerSolverBase::add_constraints() {
    int val = 0;

    return val;
}

bool PricerSolverBase::evaluate_mip_model() {
    int opt_status = model.get(GRB_IntAttr_Status);

    double objval = 0;
    if (opt_status == GRB_OPTIMAL) {
        objval = model.get(GRB_DoubleAttr_ObjVal);
        update_UB(objval);
        return objval < UB;
    } else if (opt_status == GRB_INF_OR_UNBD) {
        return false;
    } else if (opt_status == GRB_INFEASIBLE) {
        return false;
    } else if (opt_status == GRB_UNBOUNDED) {
        return false;
    } else {
        return true;
    }
}

void PricerSolverBase::remove_constraints(int first, int nb_del) {
    reformulation_model.delete_constraints(first, nb_del);
}

void PricerSolverBase::update_rows_coeff([[maybe_unused]] size_t first) {}

size_t PricerSolverBase::print_num_paths() {
    return 0UL;
}

double PricerSolverBase::get_UB() {
    return UB;
}

void PricerSolverBase::update_UB(double _ub) {
    if (_ub < UB) {
        UB = _ub;
    }
}

int PricerSolverBase::get_int_attr_model(enum MIP_Attr c) {
    int val = -1;
    switch (c) {
        case MIP_Attr_Nb_Vars:
            val = model.get(GRB_IntAttr_NumVars);
            break;
        case MIP_Attr_Nb_Constr:
            val = model.get(GRB_IntAttr_NumConstrs);
            break;
        case MIP_Attr_Status:
            val = model.get(GRB_IntAttr_Status);
            break;
        default:
            break;
    }

    return val;
}

double PricerSolverBase::get_dbl_attr_model(enum MIP_Attr c) {
    double val = -1.0;
    int    status = model.get(GRB_IntAttr_Status);
    if (status != GRB_INF_OR_UNBD && status != GRB_INFEASIBLE &&
        status != GRB_UNBOUNDED) {
        switch (c) {
            case MIP_Attr_Obj_Bound:
                val = model.get(GRB_DoubleAttr_ObjBound);
                break;
            case MIP_Attr_Obj_Bound_LP:
                val = model.get(GRB_DoubleAttr_ObjBoundC);
                break;
            case MIP_Attr_Mip_Gap:
                val = model.get(GRB_DoubleAttr_MIPGap);
                break;
            case MIP_Attr_Run_Time:
                val = model.get(GRB_DoubleAttr_Runtime);
                break;
            case MIP_Attr_Nb_Simplex_Iter:
                val = model.get(GRB_DoubleAttr_IterCount);
                break;
            case MIP_Attr_Nb_Nodes:
                val = model.get(GRB_DoubleAttr_NodeCount);
                break;
            default:
                val = std::numeric_limits<double>::max();
                break;
        }
    } else {
        switch (c) {
            case MIP_Attr_Run_Time:
                val = model.get(GRB_DoubleAttr_Runtime);
                break;
            case MIP_Attr_Nb_Simplex_Iter:
                val = model.get(GRB_DoubleAttr_IterCount);
                break;
            case MIP_Attr_Nb_Nodes:
                val = model.get(GRB_DoubleAttr_NodeCount);
                break;
            default:
                val = std::numeric_limits<double>::max();
                break;
        }
    }
    return val;
}

double PricerSolverBase::compute_reduced_cost(const OptimalSolution<>& sol,
                                              double*                  pi,
                                              double*                  lhs) {
    double result = sol.cost;
    // auto      nb_constraints = reformulation_model.get_nb_constraints();
    std::span aux_lhs{lhs, reformulation_model.size()};
    std::span aux_pi{pi, reformulation_model.size()};
    std::ranges::fill(aux_lhs, 0.0);

    for (auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        for (const auto&& [c, constr] :
             reformulation_model | ranges::views::enumerate) {
            if (c == convex_constr_id) {
                continue;
            }
            auto coeff = (*constr)(k);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * aux_pi[c];
                aux_lhs[c] += coeff;
            }
        }
    }

    double          dual = aux_pi[convex_constr_id];
    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    double          coeff = (*constr)(k);
    result -= coeff * dual;
    aux_lhs[convex_constr_id] += coeff;

    return result;
}

double PricerSolverBase::compute_lagrange(const OptimalSolution<>&   sol,
                                          const std::vector<double>& pi) {
    double result = sol.cost;
    double dual_bound = 0.0;
    // std::span aux_pi{pi, reformulation_model.size()};

    for (auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto            dual = pi[it->job];
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (fabs(coeff) > EPS_SOLVER) {
            result -= coeff * dual;
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            double dual_ = pi[c];
            double coeff_ = (*reformulation_model[c])(k);

            if (fabs(coeff_) > EPS_SOLVER) {
                result -= coeff_ * dual_;
            }
        }
    }

    result = std::min(0.0, result);

    for (const auto&& [c, constr] :
         reformulation_model | ranges::views::enumerate) {
        if (c == convex_constr_id) {
            continue;
        }

        dual_bound += constr->get_rhs() * pi[c];
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

double PricerSolverBase::compute_subgradient(const OptimalSolution<>& sol,
                                             double* subgradient) {
    std::span aux_subgradient{subgradient, reformulation_model.size()};
    auto      rhs = -reformulation_model[convex_constr_id]->get_rhs();

    for (size_t i = 0; const auto& constr : reformulation_model) {
        aux_subgradient[i] = constr->get_rhs();
        ++i;
    }

    for (auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (fabs(coeff) > EPS_SOLVER) {
            aux_subgradient[k.get_j()] -= coeff * rhs;
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            auto coeff_ = (*reformulation_model[c])(k);

            if (fabs(coeff_) > EPS_SOLVER) {
                aux_subgradient[c] -= coeff_ * rhs;
            }
        }
    }

    aux_subgradient[convex_constr_id] += rhs;

    return 0.0;
}

void PricerSolverBase::calculate_constLB(double* pi) {
    constLB = 0.0;
    std::span aux_pi{pi, reformulation_model.size()};
    for (const auto&& [c, constr] :
         reformulation_model | ranges::views::enumerate) {
        if (c == convex_constr_id) {
            continue;
        }
        constLB += constr->get_rhs() * aux_pi[c];
    }
}
/*-------------------------------------------------------------------------*/
/*------------------initialize static members------------------------------*/
/*-------------------------------------------------------------------------*/

const std::shared_ptr<GRBEnv> PricerSolverBase::genv =
    std::make_shared<GRBEnv>();
