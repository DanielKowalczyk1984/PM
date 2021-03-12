#include "PricerSolverBase.hpp"
#include <bits/ranges_algobase.h>
#include <fmt/core.h>
#include <algorithm>
#include <limits>
#include <memory>
#include "Instance.h"
#include "gurobi_c.h"

/**
 * PricerSolverBase default COnstructor
 **/
// PricerSolverBase::PricerSolverBase(GPtrArray*  _jobs,
//                                    int         _num_machines,
//                                    const char* _p_name,
//                                    double      _ub)
//     : jobs(_jobs->pdata, _jobs->len),
//       convex_constr_id(_jobs->len),
//       convex_rhs(_num_machines),
//       problem_name(_p_name),
//       env(std::make_shared<GRBEnv>()),
//       model(GRBModel(*env)),
//       reformulation_model(_jobs->len, _num_machines),
//       is_integer_solution(false),
//       constLB(0.0),
//       UB(_ub)

// {
//     try {
//         model.set(GRB_IntParam_Method, GRB_METHOD_AUTO);
//         model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//         model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
//     } catch (const GRBException& e) {
//         fmt::print("Error code = {}\n", e.getErrorCode());
//         fmt::print(e.getMessage());
//     } catch (...) {
//         fmt::print("Exception during optimization\n");
//     }
// }

PricerSolverBase::PricerSolverBase(const Instance& instance)
    : jobs(instance.jobs),
      convex_constr_id(instance.nb_jobs),
      convex_rhs(instance.nb_machines),
      problem_name(),
      env(std::make_shared<GRBEnv>()),
      model(GRBModel(*env)),
      reformulation_model(instance.nb_jobs, instance.nb_machines),
      is_integer_solution(false),
      constLB(0.0),
      UB(std::numeric_limits<int>::max()) {
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
      model(other.model),
      reformulation_model(other.reformulation_model),
      is_integer_solution(other.is_integer_solution),
      constLB(other.constLB),
      UB(other.UB) {
    model.set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
};

PricerSolverBase::~PricerSolverBase() = default;

int PricerSolverBase::add_constraints() {
    int val = 0;

    return val;
}

void PricerSolverBase::remove_constraints(int first, int nb_del) {
    reformulation_model.delete_constraints(first, nb_del);
}

void PricerSolverBase::update_rows_coeff(int first) {}

void PricerSolverBase::print_num_paths() {}

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
        for (int c = 0; const auto& constr : reformulation_model) {
            if (c == convex_constr_id) {
                continue;
            }
            auto coeff = (*constr)(k);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * aux_pi[c];
                aux_lhs[c] += coeff;
            }
            ++c;
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

double PricerSolverBase::compute_lagrange(const OptimalSolution<>& sol,
                                          double*                  pi) {
    double    result = sol.cost;
    double    dual_bound = 0.0;
    std::span aux_pi{pi, reformulation_model.size()};

    for (auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto            dual = aux_pi[it->job];
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (fabs(coeff) > EPS_SOLVER) {
            result -= coeff * dual;
        }

        for (int c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            auto*  constr_ = reformulation_model[c].get();
            double dual_ = aux_pi[c];
            double coeff_ = (*constr)(k);

            if (fabs(coeff_) > EPS_SOLVER) {
                result -= coeff_ * dual_;
            }
        }
    }

    result = std::min(0.0, result);

    for (int c = 0; const auto& constr : reformulation_model) {
        if (c == convex_constr_id) {
            continue;
        }

        dual_bound += constr->get_rhs() * aux_pi[c];
        ++c;
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

double PricerSolverBase::compute_subgradient(const OptimalSolution<>& sol,
                                             double* subgradient) {
    std::span aux_subgradient{subgradient, reformulation_model.size()};
    auto      convex_rhs = -reformulation_model[convex_constr_id]->get_rhs();

    for (size_t i = 0; const auto& constr : reformulation_model) {
        aux_subgradient[i] = constr->get_rhs();
        ++i;
    }

    for (auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (fabs(coeff) > EPS_SOLVER) {
            aux_subgradient[k.get_j()] -= coeff * convex_rhs;
        }

        for (int c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            auto* constr_ = reformulation_model[c].get();
            auto  coeff_ = (*constr)(k);

            if (fabs(coeff_) > EPS_SOLVER) {
                aux_subgradient[c] -= coeff_ * convex_rhs;
            }
        }
    }

    aux_subgradient[convex_constr_id] += convex_rhs;

    return 0.0;
}

void PricerSolverBase::calculate_constLB(double* pi) {
    constLB = 0.0;
    std::span aux_pi{pi, reformulation_model.size()};
    for (int i = 0; const auto& constr : reformulation_model) {
        if (i == convex_constr_id) {
            continue;
        }
        constLB += constr->get_rhs() * aux_pi[i];
        ++i;
    }
}

// extern "C" {
// double call_get_UB(PricerSolverBase* solver) {
//     return solver->get_UB();
// }

// void call_update_UB(PricerSolverBase* solver, double _ub) {
//     solver->update_UB(_ub);
// }

// void call_evaluate_nodes(PricerSolverBase* solver, double* pi) {
//     solver->calculate_constLB(pi);
//     solver->evaluate_nodes(pi);
// }

// int call_is_integer_solution(PricerSolverBase* solver) {
//     return solver->get_is_integer_solution();
// }
// }

// inline std::vector<BddCoeff>& PricerSolverBase::get_lp_sol() {
//     return lp_sol;
// };
