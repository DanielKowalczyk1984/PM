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

typedef struct wctlp_warmstart {
    int     row_count;
    int     column_count;
    int*    rstat;
    int*    column_status;
    double* d_norm;
} wctlp_warmstart;

#define wctlp_CONT 0
#define wctlp_BIN 1
#define wctlp_INT 2

#define wctlp_EQUAL 'E'
#define wctlp_LESS_EQUAL 'L'
#define wctlp_GREATER_EQUAL 'G'

#define wctlp_LOWER -1
#define wctlp_BASIC 0
#define wctlp_UPPER -2
#define wctlp_FREE -3

#define WCTLP_LOADED 1
#define WCTLP_OPTIMAL 2
#define WCTLP_INFEASIBLE 3
#define WCTLP_INF_OR_UNBOUNDED 4
#define WCTLP_UNBOUNDED 5
#define WCTLP_CUTOFF 6
#define WCTLP_ITERATION_LIMIT 7
#define WCTLP_NODE_LIMIT 8
#define WCTLP_TIME_LIMIT 9
#define WCTLP_SOLUTION_LIMIT 10
#define WCTLP_INTERRUPTED 11
#define WCTLP_NUMERIC 12
#define WCTLP_SUBOPTIMAL 13
#define WCTLP_IN_PROGRESS 14
#define WCTLP_USER_OBJ_LIMIT 15

#define wctlp_MIN 1
#define wctlp_MAX -1

int  wctlp_init(wctlp** lp, const char* name);
void wctlp_free(wctlp** lp);

int wctlp_optimize(wctlp* lp, int* status);
int wctlp_objval(wctlp*, double* obj);
int wctlp_pi(wctlp*, double* pi);
int wctlp_chg_upperbounds(wctlp* lp, int first, int last, double ub);
int wctlp_x(wctlp*, double* x, int first);
int wctlp_rc(wctlp*, double* rc, int first);

int wctlp_basis_cols(wctlp* lp, int* column_status, int first);
int wctlp_change_obj(wctlp* lp, int start, int len, double* values);
int wctlp_addrow(wctlp* lp, int nb_non_zero, int* column_indices, double* cval,
                 char sense, double rhs, char* name);
int wctlp_addcol(wctlp* lp, int nb_non_zero, int* column_indices, double* cval,
                 double obj, double lb, double ub, char vartype, char* name);
int wctlp_addrows(wctlp* lp, int nb_rows, int nb_zero, int* start,
                  int* column_indices, double* coeff_val, char* sense,
                  double* rhs, char** name);
int wctlp_deletecols(wctlp* lp, int first_column_ind, int last_column_ind);
int wctlp_addcols(wctlp* lp, int num_vars, int nb_zero, int* start,
                  int* row_indices, double* coeff_val, double* obj, double* lb,
                  double* ub, char* vtype, char** name);

int wctlp_set_coltypes(wctlp* lp, char sense);
int wctlp_setbound(wctlp* lp, int col, char lower_or_upper, double bound);
int wctlp_obj_sense(wctlp* lp, int sense);
int wctlp_setnodelimit(wctlp* lp, int mip_node_limit);
int wctlp_set_cutoff(wctlp* lp, double cutoff);

int  wctlp_write(wctlp* lp, const char* fname);
void wctlp_warmstart_free(wctlp_warmstart** w);
void wctlp_printerrorcode(int c);

int wctlp_status(wctlp* lp, int* status);
int wctlp_chg_lb_var(wctlp* lp, int var, double lb);
int wctlp_pi_inf(wctlp* lp, double* pi);
int wctlp_get_nb_rows(wctlp* lp, int* nb_rows);
int wctlp_get_nb_cols(wctlp* lp, int* nb_cols);
int wctlp_chgcoeff(wctlp* lp, int cnt, int* column_indices, int* var_indices,
                   double* cval);
int wctlp_getcoeff(wctlp* lp, int* column_indices, int* var_indices,
                   double* cval);

double lp_int_tolerance(void);

int wctlp_get_rhs(wctlp* lp, double* rhs);
int wctlp_compute_IIS(wctlp* lp);

#ifdef __cplusplus
}
#endif
#endif
