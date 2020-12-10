#include "job.h"
#include "lp.h"
#include "scheduleset.h"
#include "util.h"
#include "wct.h"
#include "wctprivate.h"

int grab_integer_solution(NodeData* pd, double* x, double tolerance) {
    int          val = 0;
    double       test_incumbent = .0;
    double       incumbent;
    int          tot_weighted = 0;
    ScheduleSet* tmp_schedule;
    Job*         tmp_j;

    val = lp_interface_objval(pd->RMP, &incumbent);
    CCcheck_val_2(val, "lp_interface_objval failed");
    val = lp_interface_get_nb_cols(pd->RMP, &pd->nb_cols);
    CCcheck_val_2(val, "Failed get nb_cols");

    schedulesets_free(&(pd->bestcolors), &(pd->nb_best));
    pd->bestcolors = CC_SAFE_MALLOC(pd->nb_machines, ScheduleSet);
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    pd->nb_best = 0;

    assert(pd->nb_cols == pd->localColPool->len);
    for (guint i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        test_incumbent += x[i];

        if (x[i] >= 1.0 - tolerance) {
            int j = pd->nb_best;
            scheduleset_init(pd->bestcolors + j);

            g_ptr_array_set_size(pd->bestcolors[j].job_list,
                                 tmp_schedule->job_list->len);
            for (guint k = 0; k < tmp_schedule->job_list->len; ++k) {
                tmp_j = (Job*)g_ptr_array_index(tmp_schedule->job_list, k);
                g_ptr_array_add(pd->bestcolors[j].job_list, tmp_j);
                pd->bestcolors[j].total_processing_time +=
                    tmp_j->processing_time;
                pd->bestcolors[j].total_weighted_completion_time +=
                    value_Fj(pd->bestcolors[j].total_processing_time, tmp_j);
            }

            pd->nb_best++;
            tot_weighted += pd->bestcolors[j].total_weighted_completion_time;

            if (pd->nb_best > pd->nb_machines) {
                printf(
                    "ERROR: \"Integral\" solution turned out to be not "
                    "integral!\n");
                fflush(stdout);
                val = 1;
                goto CLEAN;
            }
        }
    }

    /** Write a check function */
    printf("Intermediate schedule:\n");
    print_schedule(pd->bestcolors, pd->nb_best);
    printf("with total weight %d\n", tot_weighted);
    assert(fabs((double)tot_weighted - incumbent) <= 0.00001);

    if (tot_weighted < pd->upper_bound) {
        pd->upper_bound = tot_weighted;
        pd->best_objective = tot_weighted;
    }

    if (pd->upper_bound == pd->lower_bound) {
        pd->status = finished;
    }

CLEAN:
    return val;
}

int add_artificial_var_to_rmp(NodeData* pd) {
    int val = 0;
    int nb_jobs = pd->nb_jobs;

    for (int ii = 0; ii < nb_jobs; ii++) {
        ScheduleSet* set = scheduleset_alloc(nb_jobs);
        Job*         tmp_j = (Job*)g_ptr_array_index(pd->jobarray, ii);
        g_ptr_array_add(set->job_list, tmp_j);
        set->total_processing_time = tmp_j->processing_time;
        set->total_weighted_completion_time =
            value_Fj(tmp_j->processing_time, tmp_j);
        g_ptr_array_add(pd->localColPool, set);
    }

    return val;
}

int add_lhs_scheduleset_to_rmp(ScheduleSet* set, NodeData* pd) {
    int     val = 0;
    gsize   aux;
    double* lhs_coeff = &g_array_index(pd->lhs_coeff, double, 0);

    int* aux_int_array = g_array_steal(pd->id_row, &aux);
    CC_IFFREE(aux_int_array, int);
    double* aux_double_array = g_array_steal(pd->coeff_row, &aux);
    CC_IFFREE(aux_double_array, double);

    for (int j = 0; j < pd->nb_rows; j++) {
        if (fabs(lhs_coeff[j]) > 1e-10) {
            g_array_append_val(pd->id_row, j);
            g_array_append_val(pd->coeff_row, lhs_coeff[j]);
        }
    }

    int*    id_constraint = &g_array_index(pd->id_row, int, 0);
    double* coeff_constraint = &g_array_index(pd->coeff_row, double, 0);
    int     len = pd->id_row->len;

    val = lp_interface_get_nb_cols(pd->RMP, &(pd->nb_cols));
    CCcheck_val_2(val, "Failed to get the number of cols");
    set->id = pd->nb_cols - pd->id_pseudo_schedules;
    val = lp_interface_addcol(pd->RMP, len, id_constraint, coeff_constraint,
                              (double)set->total_weighted_completion_time, 0.0,
                              GRB_INFINITY, lp_interface_CONT, NULL);
    CCcheck_val_2(val, "Failed to add column to lp");

CLEAN:
    return val;
}

int add_scheduleset_to_rmp(ScheduleSet* set, NodeData* pd) {
    int        val = 0;
    int        row_ind;
    int        var_ind;
    double     cval;
    int        nb_jobs = pd->nb_jobs;
    GPtrArray* members = set->job_list;
    wctlp*     lp = pd->RMP;
    Job*       job;

    val = lp_interface_get_nb_cols(lp, &(pd->nb_cols));
    CCcheck_val_2(val, "Failed to get the number of cols");
    var_ind = pd->nb_cols;
    set->id = pd->nb_cols - pd->id_pseudo_schedules;
    val = lp_interface_addcol(lp, 0, NULL, NULL,
                              (double)set->total_weighted_completion_time, 0.0,
                              GRB_INFINITY, lp_interface_CONT, NULL);
    CCcheck_val_2(val, "Failed to add column to lp");

    for (unsigned i = 0; i < members->len; ++i) {
        job = (Job*)g_ptr_array_index(members, i);
        row_ind = job->job;
        val = lp_interface_getcoeff(lp, &row_ind, &var_ind, &cval);
        CCcheck_val_2(val, "Failed lp_interface_getcoeff");
        cval += 1.0;
        set->num[job->job] += 1;
        val = lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);
        CCcheck_val_2(val, "Failed lp_interface_chgcoeff");
    }

    row_ind = nb_jobs;
    cval = -1.0;
    val = lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);
    CCcheck_val_2(val, "Failed lp_interface_chgcoeff");

CLEAN:
    return val;
}

void g_add_col_to_lp(gpointer data, gpointer user_data) {
    ScheduleSet* tmp = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;
    add_scheduleset_to_rmp(tmp, pd);
}

int build_rmp(NodeData* pd) {
    int     val = 0;
    int     nb_jobs = pd->nb_jobs;
    int     nb_machines = pd->nb_machines;
    Parms*  parms = pd->parms;
    double* lb = NULL;
    double* ub = NULL;
    double* obj = NULL;
    char*   vtype = NULL;
    int*    start_vars = NULL;
    int*    rows_ind = NULL;
    double* coeff_vals = NULL;
    int*    int_values = NULL;
    double* dbl_values = NULL;
    int*    start = CC_SAFE_MALLOC(nb_jobs, int);
    double* rhs = CC_SAFE_MALLOC(nb_jobs, double);
    char*   sense = CC_SAFE_MALLOC(nb_jobs, char);
    val = lp_interface_init(&(pd->RMP), NULL);
    CCcheck_val_2(val, "lp_interface_init failed");

    fill_int(start, nb_jobs, 0);
    fill_dbl(rhs, nb_jobs, 1.0);
    fill_char(sense, nb_jobs, GRB_GREATER_EQUAL);

    /**
     * add assignment constraints
     */
    lp_interface_get_nb_rows(pd->RMP, &(pd->id_assignment_constraint));
    lp_interface_addrows(pd->RMP, nb_jobs, 0, start, NULL, NULL, sense, rhs,
                         NULL);

    /**
     * add number of machines constraint (convexification)
     */
    lp_interface_get_nb_rows(pd->RMP, &(pd->id_convex_constraint));
    val = lp_interface_addrow(pd->RMP, 0, (int*)NULL, (double*)NULL,
                              lp_interface_GREATER_EQUAL, -(double)nb_machines,
                              (char*)NULL);
    CCcheck_val_2(val, "Failed to add convexification constraint");
    lp_interface_get_nb_rows(pd->RMP, &(pd->id_valid_cuts));
    pd->nb_rows = pd->id_valid_cuts;

    /**
     * construct artificial variables in RMP
     */
    lp_interface_get_nb_cols(pd->RMP, &(pd->id_art_var_assignment));
    pd->id_art_var_convex = nb_jobs;
    pd->id_art_var_cuts = nb_jobs + 1;
    pd->id_next_var_cuts = pd->id_art_var_cuts;
    int nb_vars = nb_jobs + 1 + pd->max_nb_cuts;
    pd->id_pseudo_schedules = nb_vars;

    lb = CC_SAFE_MALLOC(nb_vars, double);
    CCcheck_NULL_2(lb, "Failed to allocate memory");
    ub = CC_SAFE_MALLOC(nb_vars, double);
    CCcheck_NULL_2(ub, "Failed to allocate memory");
    obj = CC_SAFE_MALLOC(nb_vars, double);
    CCcheck_NULL_2(obj, "Failed to allocate memory");
    vtype = CC_SAFE_MALLOC(nb_vars, char);
    CCcheck_NULL_2(vtype, "Failed to allocate memory");
    start_vars = CC_SAFE_MALLOC(nb_vars + 1, int);
    CCcheck_NULL_2(start_vars, "Failed to allocate memory");

    for (int j = 0; j < nb_vars; j++) {
        lb[j] = 0.0;
        ub[j] = GRB_INFINITY;
        vtype[j] = GRB_CONTINUOUS;
        start_vars[j] = nb_jobs + 1;
        obj[j] = (pd->opt_sol->tw + 1) * 100;
    }
    start_vars[nb_vars] = nb_jobs + 1;

    rows_ind = CC_SAFE_MALLOC(nb_jobs + 1, int);
    CCcheck_NULL_2(rows_ind, "Failed to allocate memory");
    coeff_vals = CC_SAFE_MALLOC(nb_jobs + 1, double);
    CCcheck_NULL_2(coeff_vals, "Failed to allocate memory");

    for (int j = 0; j < nb_jobs + 1; j++) {
        rows_ind[j] = j;
        start_vars[j] = j;
        coeff_vals[j] = 1.0;
    }

    val = lp_interface_addcols(pd->RMP, nb_vars, nb_jobs + 1, start_vars,
                               rows_ind, coeff_vals, obj, lb, ub, vtype, NULL);
    CCcheck_val(val, "failed construct cols");
    val = lp_interface_get_nb_cols(pd->RMP, &(pd->id_pseudo_schedules));

    /** add columns from localColPool */
    // add_artificial_var_to_rmp(pd);
    g_ptr_array_foreach(pd->localColPool, g_add_col_to_lp, pd);

    /**
     * Some aux variables for column generation
     */

    dbl_values = CC_SAFE_MALLOC(pd->nb_rows, double);
    CCcheck_NULL_2(dbl_values, "Failed to allocate memory");
    fill_dbl(dbl_values, pd->nb_rows, 0.0);
    int_values = CC_SAFE_MALLOC(pd->nb_rows, int);
    CCcheck_NULL_2(int_values, "Failed to allocate memory");
    fill_int(int_values, pd->nb_rows, 0);

    pd->pi = g_array_sized_new(FALSE, FALSE, sizeof(double), pd->nb_rows);
    CCcheck_NULL_2(pd->pi, "Failed to allocate memory");
    g_array_append_vals(pd->pi, dbl_values, pd->nb_rows);
    pd->slack = g_array_sized_new(FALSE, FALSE, sizeof(double), pd->nb_rows);
    CCcheck_NULL_2(pd->slack, "Failed to allocate memory");
    g_array_append_vals(pd->slack, dbl_values, pd->nb_rows);
    pd->rhs = g_array_sized_new(FALSE, FALSE, sizeof(double), pd->nb_rows);
    CCcheck_NULL_2(pd->rhs, "Failed to allocate memory");
    g_array_append_vals(pd->rhs, dbl_values, pd->nb_rows);
    val = lp_interface_get_rhs(pd->RMP, &g_array_index(pd->rhs, double, 0));
    CCcheck_val_2(val, "Failed to get RHS");
    pd->lhs_coeff =
        g_array_sized_new(FALSE, FALSE, sizeof(double), pd->nb_rows);
    CCcheck_NULL_2(pd->lhs_coeff, "Failed to allocate memory");
    g_array_append_vals(pd->lhs_coeff, dbl_values, pd->nb_rows);
    pd->id_row = g_array_sized_new(FALSE, FALSE, sizeof(int), pd->nb_rows);
    CCcheck_NULL_2(pd->id_row, "Failed to allocate memory");
    g_array_append_vals(pd->id_row, int_values, pd->nb_rows);
    pd->coeff_row =
        g_array_sized_new(FALSE, FALSE, sizeof(double), pd->nb_rows);
    CCcheck_NULL_2(pd->coeff_row, "Failed to allocate memory");
    g_array_append_vals(pd->coeff_row, dbl_values, pd->nb_rows);

    pd->eta_in = 0.0;
    pd->eta_out = 0.0;
    if (parms->pricing_solver < dp_solver) {
        pd->x_e = CC_SAFE_MALLOC(2 * get_nb_vertices(pd->solver), double);
        CCcheck_NULL_2(pd->x_e, "Failed to allocate memory");
    }

CLEAN:

    if (val) {
        lp_interface_free(&(pd->RMP));
        g_array_free(pd->pi, TRUE);
        // g_array_free(pd->pi_in, TRUE);
        // g_array_free(pd->pi_out, TRUE);
        // g_array_free(pd->pi_sep, TRUE);
        // g_array_free(pd->subgradient, TRUE);
        // g_array_free(pd->subgradient_in, TRUE);
        g_array_free(pd->rhs, TRUE);
        g_array_free(pd->lhs_coeff, TRUE);
    }
    CC_IFFREE(sense, char);
    CC_IFFREE(rhs, double);
    CC_IFFREE(start, int);
    CC_IFFREE(lb, double);
    CC_IFFREE(ub, double);
    CC_IFFREE(vtype, char);
    CC_IFFREE(obj, double);
    CC_IFFREE(start_vars, int);
    CC_IFFREE(rows_ind, int);
    CC_IFFREE(coeff_vals, double);
    CC_IFFREE(dbl_values, double);
    CC_IFFREE(int_values, int);
    return val;
}

int get_solution_lp_lowerbound(NodeData* pd) {
    int  val = 0;
    Job* tmp_j;

    val = lp_interface_get_nb_cols(pd->RMP, &pd->nb_cols);
    CCcheck_val_2(val, "Failed in wctlp_get_nb_cols");
    pd->lambda = CC_SAFE_REALLOC(pd->lambda, pd->nb_cols, double);
    CCcheck_NULL_2(pd->lambda, "Failed to allocate memory");
    assert(pd->nb_cols == pd->localColPool->len);
    lp_interface_x(pd->RMP, pd->lambda, 0);

    for (int i = 0; i < pd->nb_cols; ++i) {
        ScheduleSet* tmp =
            ((ScheduleSet*)g_ptr_array_index(pd->localColPool, i));
        if (pd->lambda[i]) {
            printf("%f: ", pd->lambda[i]);
            g_ptr_array_foreach(tmp->job_list, g_print_machine, NULL);
            printf("\n");
            for (unsigned j = 0; j < tmp->job_list->len; ++j) {
                tmp_j = (Job*)g_ptr_array_index(tmp->job_list, j);
                printf("%d (%d) ", tmp->num[tmp_j->job],
                       tmp_j->processing_time);
            }
            printf("\n");
        }
    }

CLEAN:
    return val;
}
