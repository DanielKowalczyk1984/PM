#include "PricerSolverBase.hpp"
#include <fmt/core.h>  // for print
#include <algorithm>   // for min, __fill_fn
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>                                       // for fabs
#include <cstddef>                                     // for size_t
#include <ext/alloc_traits.h>                          // for __alloc_traits...
#include <limits>                                      // for numeric_limits
#include <memory>                                      // for __shared_ptr_a...
#include <range/v3/iterator/basic_iterator.hpp>        // for basic_iterator
#include <range/v3/iterator/unreachable_sentinel.hpp>  // for operator==
#include <range/v3/range/conversion.hpp>               // for to_container::fn
#include <range/v3/view/enumerate.hpp>                 // for enumerate_fn
#include <range/v3/view/transform.hpp>                 // for trqnsform_view
#include <range/v3/view/view.hpp>                      // for operator|
#include <range/v3/view/zip.hpp>                       // for zip_view
#include <range/v3/view/zip_with.hpp>                  // for iter_zip_with_...
#include <span>                                        // for span
#include <vector>                                      // for vector
#include "Column.h"                                    // for Column
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
        fmt::print("{}", e.getMessage());
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
    switch (opt_status) {
        case GRB_OPTIMAL:
            objval = model.get(GRB_DoubleAttr_ObjVal);
            update_UB(objval);
            return objval < UB;
            break;
        case GRB_INFEASIBLE:
        case GRB_INF_OR_UNBD:
        case GRB_UNBOUNDED:
            return false;
            break;
        default:
            return true;
    }
}

bool PricerSolverBase::compute_sub_optimal_duals(
    const double*                               lambda,
    const std::vector<std::shared_ptr<Column>>& columns) {
    GRBModel sub_optimal(*genv);
    // sub_optimal.set(GRB_IntParam_OutputFlag, 1);
    sub_optimal.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    std::vector<GRBVar> beta{jobs.size()};
    std::vector<GRBVar> eta;
    double              LB{};
    auto                removed = false;

    for (auto& it : beta) {
        it = sub_optimal.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
    }

    auto last = sub_optimal.addVar(
        0.0, GRB_INFINITY, -static_cast<double>(convex_rhs), GRB_CONTINUOUS);

    std::span<const double> aux_cols{lambda, columns.size()};

    for (auto&& [set, x] : ranges::views::zip(columns, aux_cols)) {
        if (x > EPS_SOLVER) {
            // auto* tmp = columns[i].get();
            eta.emplace_back(sub_optimal.addVar(0.0, GRB_INFINITY, 1.0, 'C'));
            GRBLinExpr expr = -last;
            for (auto& it : set->job_list) {
                expr += beta[it->job];
            }
            expr += eta.back();
            sub_optimal.addConstr(expr, '=',
                                  set->total_weighted_completion_time);
            LB += set->total_weighted_completion_time * x;
        } else {
            // auto*      tmp = columns[i].get();
            GRBLinExpr expr = -last;
            for (auto& it : set->job_list) {
                expr += beta[it->job];
            }
            sub_optimal.addConstr(expr, '<',
                                  set->total_weighted_completion_time);
        }
    }

    GRBLinExpr expr = -static_cast<double>(convex_rhs) * last;
    for (auto& it : beta) {
        expr += it;
    }
    sub_optimal.addConstr(expr, '>', LB - RC_FIXING);

    auto cont = false;
    do {
        sub_optimal.update();
        sub_optimal.optimize();

        auto pi = beta |
                  ranges::views::transform([&](const auto& tmp) -> double {
                      return tmp.get(GRB_DoubleAttr_X);
                  }) |
                  ranges::to_vector;
        pi.push_back(last.get(GRB_DoubleAttr_X));

        auto sol = pricing_algorithm(pi.data());
        auto rc = sol.cost + pi.back();
        for (auto& it : sol.jobs) {
            rc -= pi[it->job];
        }

        if (rc < -RC_FIXING) {
            GRBLinExpr expr_pricing = -last;
            for (auto& it : sol.jobs) {
                expr_pricing += beta[it->job];
            }
            sub_optimal.addConstr(expr_pricing, '<', sol.cost);
            cont = true;
        } else {
            cont = false;
            calculate_constLB(pi.data());
            removed = evaluate_nodes(pi.data());
        }

    } while (cont);

    return removed;
}

void PricerSolverBase::remove_constraints(int first, int nb_del) {
    reformulation_model.delete_constraints(first, nb_del);
}

void PricerSolverBase::update_rows_coeff([[maybe_unused]] size_t first) {}

boost::multiprecision::cpp_int PricerSolverBase::print_num_paths() {
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

double PricerSolverBase::compute_reduced_cost(const PricingSolution<>& sol,
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

double PricerSolverBase::compute_reduced_cost_simple(
    const PricingSolution<>& sol,
    double*                  pi) {
    double result = sol.cost;
    // auto      nb_constraints = reformulation_model.get_nb_constraints();
    std::span aux_pi{pi, reformulation_model.size()};

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
            }
        }
    }

    double          dual = aux_pi[convex_constr_id];
    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    double          coeff = (*constr)(k);
    result -= coeff * dual;

    return result;
}

double PricerSolverBase::compute_lagrange(const PricingSolution<>&   sol,
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

double PricerSolverBase::compute_subgradient(const PricingSolution<>& sol,
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
/*------------------initialize static
 * members------------------------------*/
/*-------------------------------------------------------------------------*/

const std::shared_ptr<GRBEnv> PricerSolverBase::genv =
    std::make_shared<GRBEnv>();
