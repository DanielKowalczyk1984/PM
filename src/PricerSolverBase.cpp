#include "PricerSolverBase.hpp"

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
      model(new GRBModel(*env)) { }

PricerSolverBase::PricerSolverBase(GPtrArray* _jobs, int _num_machines,
                                   GPtrArray* _ordered_jobs, const char* p_name)
    : jobs(_jobs),
      nb_jobs(_jobs->len),
      num_machines(_num_machines),
      ordered_jobs(_ordered_jobs),
      nb_layers(ordered_jobs->len),
      problem_name(p_name),
      env(new GRBEnv()),
      model(new GRBModel(*env)) {}

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
            break;
    }

    return val;
}
