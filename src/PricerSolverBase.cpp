#include "PricerSolverBase.hpp"
#include <limits>

/**
 * PricerSolverBase default COnstructor
 **/
PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   const char* p_name)
    : jobs(_jobs),
      nb_jobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(nullptr),
      nb_layers(0),
      problem_name(p_name),
      env(new GRBEnv()),
      model(new GRBModel(*env)),
      reformulation_model(jobs->len, _num_machines)
       {
    model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->set(GRB_IntParam_Presolve, 2);
    model->set(GRB_DoubleParam_MIPGap, 1e-6);
    model->set(GRB_DoubleParam_TimeLimit, 1800);
}

PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   GPtrArray* _ordered_jobs, const char* p_name)
    : jobs(_jobs),
      nb_jobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(_ordered_jobs),
      nb_layers(ordered_jobs->len),
      problem_name(p_name),
      env(new GRBEnv()),
      model(new GRBModel(*env)),
      reformulation_model(jobs->len, _num_machines)
       {
    model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->set(GRB_IntParam_Presolve, 2);
    model->set(GRB_DoubleParam_MIPGap, 0.0);
    model->set(GRB_DoubleParam_TimeLimit, 1800);
}

PricerSolverBase::~PricerSolverBase() {}

void PricerSolverBase::print_num_paths() {}

int PricerSolverBase::get_int_attr_model(enum MIP_Attr c) {
    int val = -1;
    switch (c) {
        case MIP_Attr_Nb_Vars:
            val = model->get(GRB_IntAttr_NumVars);
            break;
        case MIP_Attr_Nb_Constr:
            val = model->get(GRB_IntAttr_NumConstrs);
            break;
        case MIP_Attr_Status:
            val = model->get(GRB_IntAttr_Status);
            break;
        default:
            break;
    }

    return val;
}

double PricerSolverBase::get_dbl_attr_model(enum MIP_Attr c) {
    double val = -1.0;
    int    status = model->get(GRB_IntAttr_Status);
    if (status != GRB_INF_OR_UNBD && status != GRB_INFEASIBLE &&
        status != GRB_UNBOUNDED) {
        switch (c) {
            case MIP_Attr_Obj_Bound:
                val = model->get(GRB_DoubleAttr_ObjBound);
                break;
            case MIP_Attr_Obj_Bound_LP:
                val = model->get(GRB_DoubleAttr_ObjBoundC);
                break;
            case MIP_Attr_Mip_Gap:
                val = model->get(GRB_DoubleAttr_MIPGap);
                break;
            case MIP_Attr_Run_Time:
                val = model->get(GRB_DoubleAttr_Runtime);
                break;
            case MIP_Attr_Nb_Simplex_Iter:
                val = model->get(GRB_DoubleAttr_IterCount);
                break;
            case MIP_Attr_Nb_Nodes:
                val = model->get(GRB_DoubleAttr_NodeCount);
                break;
            default:
                val = std::numeric_limits<double>::max();
                break;
        }
    } else {
        switch (c) {
            case MIP_Attr_Run_Time:
                val = model->get(GRB_DoubleAttr_Runtime);
                break;
            case MIP_Attr_Nb_Simplex_Iter:
                val = model->get(GRB_DoubleAttr_IterCount);
                break;
            case MIP_Attr_Nb_Nodes:
                val = model->get(GRB_DoubleAttr_NodeCount);
                break;
            default:
                val = std::numeric_limits<double>::max();
                break;
        }
    }
    return val;
}

double PricerSolverBase::compute_reduced_cost(const OptimalSolution<>& sol, double *pi, double *lhs){
    double result = sol.cost;
    std::fill(lhs, lhs+reformulation_model.get_nb_constraints() , 0.0);

    for (guint j = 0; j < sol.jobs->len; j++) {
        Job* tmp_j = (Job*) g_ptr_array_index(sol.jobs, j);
        VariableKeyBase k(tmp_j->job,0);
        for(int c = 0; c < nb_jobs; c++) {
            double dual = pi[c];
            ConstraintBase *constr = reformulation_model.get_constraint(c);
            double coeff = constr->get_var_coeff(&k);

            if(fabs(coeff) > 1e-10) {
                result -= coeff*dual;
                lhs[c] += coeff;
            }
        }
    }

    double dual = pi[nb_jobs];
    ConstraintBase *constr = reformulation_model.get_constraint(nb_jobs);
    VariableKeyBase k(0,0,true);
    double coeff = constr->get_var_coeff(&k);
    result -= coeff*dual;
    lhs[nb_jobs] += coeff;

    return result;
}

double PricerSolverBase::compute_lagrange(const OptimalSolution<> &sol, double *pi) {
    double result = sol.cost;
    double dual_bound = 0.0;

    for (guint j = 0; j < sol.jobs->len; j++) {
        Job* tmp_j = (Job*) g_ptr_array_index(sol.jobs, j);
        VariableKeyBase k(tmp_j->job,0);
        for(int c = 0; c < nb_jobs; c++) {
            double dual = pi[c];
            ConstraintBase *constr = reformulation_model.get_constraint(c);
            double coeff = constr->get_var_coeff(&k);

            if(fabs(coeff) > 1e-10) {
                result -= coeff*dual;
            }
        }
    }

    result = CC_MIN(0, result);

    for(int c = 0; c < nb_jobs; c++) {
        double dual = pi[c];
        ConstraintBase *constr = reformulation_model.get_constraint(c);
        double rhs = constr->get_rhs();

        dual_bound += rhs*dual;
    }

    result = -reformulation_model.get_constraint(nb_jobs)->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

