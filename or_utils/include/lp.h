////////////////////////////////////////////////////////////////
//                                                            //
//  lp.h                                                      //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef __LP_h
#define __LP_h

#ifdef __cplusplus
extern "C" {
#endif

#include <gurobi_c.h>

#define CHECK_VAL_GRB(val, msg, env)                                    \
    {                                                                   \
        if ((val)) {                                                    \
            fprintf(stderr, "%s at %s, line %d: %s\n", (msg), __FILE__, \
                    __LINE__, GRBgeterrormsg(env));                     \
        }                                                               \
    }

#define CHECK_VAL_GRB2(val, msg, env)                                   \
    {                                                                   \
        if ((val)) {                                                    \
            fprintf(stderr, "%s at %s, line %d: %s\n", (msg), __FILE__, \
                    __LINE__, GRBgeterrormsg(env));                     \
            goto CLEAN;                                                 \
        }                                                               \
    }

typedef struct lp wctlp;

struct lp {
    GRBenv*   env;
    GRBmodel* model;

    double dbl_cutoff;
};

typedef struct lp_interface_warmstart {
    int     row_count;
    int     column_count;
    int*    rstat;
    int*    column_status;
    double* d_norm;
} lp_interface_warmstart;

#define lp_interface_CONT 0
#define lp_interface_BIN 1
#define lp_interface_INT 2

#define lp_interface_EQUAL 'E'
#define lp_interface_LESS_EQUAL 'L'
#define lp_interface_GREATER_EQUAL 'G'

#define lp_interface_LOWER -1
#define lp_interface_BASIC 0
#define lp_interface_UPPER -2
#define lp_interface_FREE -3

#define LP_INTERFACE_LOADED 1
#define LP_INTERFACE_OPTIMAL 2
#define LP_INTERFACE_INFEASIBLE 3
#define LP_INTERFACE_INF_OR_UNBOUNDED 4
#define LP_INTERFACE_UNBOUNDED 5
#define LP_INTERFACE_CUTOFF 6
#define LP_INTERFACE_ITERATION_LIMIT 7
#define LP_INTERFACE_NODE_LIMIT 8
#define LP_INTERFACE_TIME_LIMIT 9
#define LP_INTERFACE_SOLUTION_LIMIT 10
#define LP_INTERFACE_INTERRUPTED 11
#define LP_INTERFACE_NUMERIC 12
#define LP_INTERFACE_SUBOPTIMAL 13
#define LP_INTERFACE_IN_PROGRESS 14
#define LP_INTERFACE_USER_OBJ_LIMIT 15

#define lp_interface_MIN 1
#define lp_interface_MAX -1

int  lp_interface_init(wctlp** lp, const char* name);
void lp_interface_free(wctlp** lp);

int lp_interface_optimize(wctlp* lp, int* status);
int lp_interface_objval(wctlp*, double* obj);
int lp_interface_pi(wctlp*, double* pi);
int lp_interface_slack(wctlp* lp, double* slack);
int lp_interface_chg_upperbounds(wctlp* lp, int first, int last, double ub);
int lp_interface_x(wctlp*, double* x, int first);
int lp_interface_rc(wctlp*, double* rc, int first);

int lp_interface_basis_cols(wctlp* lp, int* column_status, int first);
int lp_interface_change_obj(wctlp* lp, int start, int len, double* values);
int lp_interface_addrow(wctlp*  lp,
                        int     nb_non_zero,
                        int*    column_indices,
                        double* cval,
                        char    sense,
                        double  rhs,
                        char*   name);
int lp_interface_addcol(wctlp*  lp,
                        int     nb_non_zero,
                        int*    column_indices,
                        double* cval,
                        double  obj,
                        double  lb,
                        double  ub,
                        char    vartype,
                        char*   name);
int lp_interface_addrows(wctlp*  lp,
                         int     nb_rows,
                         int     nb_zero,
                         int*    start,
                         int*    column_indices,
                         double* coeff_val,
                         char*   sense,
                         double* rhs,
                         char**  name);
int lp_interface_deletecols(wctlp* lp,
                            int    first_column_ind,
                            int    last_column_ind);
int lp_interface_deleterows(wctlp* lp, int first, int last);
int lp_interface_addcols(wctlp*  lp,
                         int     num_vars,
                         int     nb_zero,
                         int*    start,
                         int*    row_indices,
                         double* coeff_val,
                         double* obj,
                         double* lb,
                         double* ub,
                         char*   vtype,
                         char**  name);

int lp_interface_set_coltypes(wctlp* lp, char sense);
int lp_interface_setbound(wctlp* lp,
                          int    col,
                          char   lower_or_upper,
                          double bound);
int lp_interface_obj_sense(wctlp* lp, int sense);
int lp_interface_setnodelimit(wctlp* lp, int mip_node_limit);
int lp_interface_set_cutoff(wctlp* lp, double cutoff);

int  lp_interface_write(wctlp* lp, const char* fname);
void lp_interface_warmstart_free(lp_interface_warmstart** w);
void lp_interface_printerrorcode(int c);

int lp_interface_status(wctlp* lp, int* status);
int lp_interface_chg_lb_var(wctlp* lp, int var, double lb);
int lp_interface_pi_inf(wctlp* lp, double* pi);
int lp_interface_get_nb_rows(wctlp* lp, int* nb_rows);
int lp_interface_get_nb_cols(wctlp* lp, int* nb_cols);
int lp_interface_chgcoeff(wctlp*  lp,
                          int     cnt,
                          int*    column_indices,
                          int*    var_indices,
                          double* cval);
int lp_interface_getcoeff(wctlp*  lp,
                          int*    column_indices,
                          int*    var_indices,
                          double* cval);

double lp_int_tolerance(void);

int lp_interface_get_rhs(wctlp* lp, double* rhs);
int lp_interface_compute_IIS(wctlp* lp);

#ifdef __cplusplus
}
#endif
#endif
