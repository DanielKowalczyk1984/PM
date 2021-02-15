#include <fmt/core.h>
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>
#include "gurobi_c.h"
#include "job.h"
#include "lp.h"
#include "scheduleset.h"
#include "util.h"
#include "wctprivate.h"

int grab_integer_solution(NodeData* pd, double* x, double tolerance) {
    int          val = 0;
    double       incumbent = 0.0;
    int          tot_weighted = 0;
    ScheduleSet* tmp_schedule = nullptr;
    Job*         tmp_j = nullptr;

    val = lp_interface_objval(pd->RMP, &incumbent);
    val = lp_interface_get_nb_cols(pd->RMP, &pd->nb_cols);
    assert(pd->nb_cols - pd->id_pseudo_schedules == pd->localColPool->len);

    g_ptr_array_free(pd->best_schedule, TRUE);
    pd->best_schedule = g_ptr_array_new_with_free_func(g_scheduleset_free);

    for (guint i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule =
            static_cast<ScheduleSet*>(g_ptr_array_index(pd->localColPool, i));

        if (x[i + pd->id_pseudo_schedules] >= 1.0 - tolerance) {
            ScheduleSet* aux = CC_SAFE_MALLOC(1, ScheduleSet);
            scheduleset_init_bis(aux);
            aux->job_list =
                g_ptr_array_copy(tmp_schedule->job_list, NULL, NULL);
            scheduleset_recalculate(aux);
            g_ptr_array_add(pd->best_schedule, aux);

            tot_weighted += aux->total_weighted_completion_time;

            if (pd->best_schedule->len > pd->nb_machines) {
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
    fmt::print("Intermediate schedule:\n");
    // print_schedule(pd->bestcolors, pd->nb_best);
    fmt::print("with total weight {}\n", tot_weighted);
    assert(fabs((double)tot_weighted - incumbent) <= EPS);

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

int NodeData::add_artificial_var_to_rmp() {
    int val = 0;

    for (int ii = 0; ii < nb_jobs; ii++) {
        ScheduleSet* set = scheduleset_alloc(nb_jobs);
        Job*         tmp_j = static_cast<Job*>(g_ptr_array_index(jobarray, ii));
        g_ptr_array_add(set->job_list, tmp_j);
        set->total_processing_time = tmp_j->processing_time;
        set->total_weighted_completion_time =
            value_Fj(tmp_j->processing_time, tmp_j);
        g_ptr_array_add(localColPool, set);
    }

    return val;
}

int NodeData::add_lhs_scheduleset_to_rmp(ScheduleSet* set) {
    int   val = 0;
    gsize aux = 0;
    // double* lhs_coeff_tmp =
    //     static_cast<double*>(&g_array_index(lhs_coeff, double, 0));

    // int* aux_int_array = static_cast<int*>(g_array_steal(id_row, &aux));
    // CC_IFFREE(aux_int_array, int);
    id_row.clear();
    coeff_row.clear();
    // double* aux_double_array =
    //     static_cast<double*>(g_array_steal(coeff_row, &aux));
    // CC_IFFREE(aux_double_array, double);

    for (int j = 0; j < nb_rows; j++) {
        if (std::fabs(lhs_coeff[j]) > EPS_BOUND) {
            id_row.emplace_back(j);
            coeff_row.emplace_back(lhs_coeff[j]);
        }
    }

    // int*    id_constraint = static_cast<int*>(&g_array_index(id_row, int,
    // 0));
    // double* coeff_constraint =
    //     static_cast<double*>(&g_array_index(coeff_row, double, 0));
    int len = id_row.size();

    val = lp_interface_get_nb_cols(RMP, &(nb_cols));
    set->id = nb_cols - id_pseudo_schedules;
    auto cost = static_cast<double>(set->total_weighted_completion_time);
    val = lp_interface_addcol(RMP, len, id_row.data(), coeff_row.data(), cost,
                              0.0, GRB_INFINITY, lp_interface_CONT, NULL);

    return val;
}

int NodeData::add_scheduleset_to_rmp(ScheduleSet* set) {
    int        val = 0;
    int        row_ind = 0;
    int        var_ind = 0;
    double     cval = 0.0;
    double     cost = static_cast<double>(set->total_weighted_completion_time);
    GPtrArray* members = set->job_list;
    wctlp*     lp = RMP;
    Job*       job = nullptr;

    val = lp_interface_get_nb_cols(lp, &(nb_cols));
    var_ind = nb_cols;
    set->id = nb_cols - id_pseudo_schedules;
    val = lp_interface_addcol(lp, 0, nullptr, nullptr, cost, 0.0, GRB_INFINITY,
                              lp_interface_CONT, nullptr);

    for (unsigned i = 0; i < members->len; ++i) {
        job = static_cast<Job*>(g_ptr_array_index(members, i));
        row_ind = job->job;
        val = lp_interface_getcoeff(lp, &row_ind, &var_ind, &cval);
        cval += 1.0;
        val = lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);
    }

    row_ind = nb_jobs;
    cval = -1.0;
    val = lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);

    return val;
}

void g_add_col_to_lp(gpointer data, gpointer user_data) {
    ScheduleSet* tmp = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;
    pd->add_scheduleset_to_rmp(tmp);
}

int NodeData::build_rmp() {
    int                 val = 0;
    double*             dbl_values{};
    int*                int_values{};
    std::vector<int>    start(nb_jobs, 0);
    std::vector<double> rhs_tmp(nb_jobs, 1.0);
    std::vector<char>   sense(nb_jobs, GRB_GREATER_EQUAL);
    val = lp_interface_init(&(RMP), NULL);

    /**
     * add assignment constraints
     */
    lp_interface_get_nb_rows(RMP, &(id_assignment_constraint));
    lp_interface_addrows(RMP, nb_jobs, 0, start.data(), NULL, NULL,
                         sense.data(), rhs_tmp.data(), NULL);

    /**
     * add number of machines constraint (convexification)
     */
    lp_interface_get_nb_rows(RMP, &(id_convex_constraint));
    val = lp_interface_addrow(RMP, 0, (int*)NULL, (double*)NULL,
                              lp_interface_GREATER_EQUAL, -(double)nb_machines,
                              (char*)NULL);
    // CCcheck_val_2(val, "Failed to add convexification constraint");
    lp_interface_get_nb_rows(RMP, &(id_valid_cuts));
    nb_rows = id_valid_cuts;

    /**
     * construct artificial variables in RMP
     */
    lp_interface_get_nb_cols(RMP, &(id_art_var_assignment));
    id_art_var_convex = nb_jobs;
    id_art_var_cuts = nb_jobs + 1;
    id_next_var_cuts = id_art_var_cuts;
    int nb_vars = nb_jobs + 1 + max_nb_cuts;
    id_pseudo_schedules = nb_vars;

    std::vector<double> lb(nb_vars, 0.0);
    std::vector<double> ub(nb_vars, GRB_INFINITY);
    std::vector<double> obj(nb_vars, 100.0 * (opt_sol->tw + 1));
    std::vector<char>   vtype(nb_vars, GRB_CONTINUOUS);
    std::vector<int>    start_vars(nb_vars + 1, nb_jobs + 1);

    start_vars[nb_vars] = nb_jobs + 1;

    std::vector<int> rows_ind(nb_jobs + 1);
    std::iota(rows_ind.begin(), rows_ind.end(), 0);
    std::vector<double> coeff_vals(nb_jobs + 1, 1.0);
    // std::iota(start_vars.begin(), start_vars.begin() + nb_jobs + 1, 0);

    for (int j = 0; j < nb_jobs + 1; j++) {
        //     rows_ind[j] = j;
        start_vars[j] = j;
        //     coeff_vals[j] = 1.0;
    }

    val = lp_interface_addcols(RMP, nb_vars, nb_jobs + 1, start_vars.data(),
                               rows_ind.data(), coeff_vals.data(), obj.data(),
                               lb.data(), ub.data(), vtype.data(), nullptr);
    // CCcheck_val(val, "failed construct cols");
    val = lp_interface_get_nb_cols(RMP, &(id_pseudo_schedules));

    /** add columns from localColPool */
    // add_artificial_var_to_rmp(pd);
    g_ptr_array_foreach(localColPool, g_add_col_to_lp, this);

    /**
     * Some aux variables for column generation
     */

    dbl_values = CC_SAFE_MALLOC(nb_rows, double);
    fill_dbl(dbl_values, nb_rows, 0.0);
    int_values = CC_SAFE_MALLOC(nb_rows, int);
    fill_int(int_values, nb_rows, 0);

    pi.resize(nb_rows, 0.0);
    // g_array_append_vals(pi, dbl_values, nb_rows);
    slack.resize(nb_rows, 0.0);
    // g_array_append_vals(slack, dbl_values, nb_rows);
    // rhs = g_array_sized_new(FALSE, FALSE, sizeof(double), nb_rows);
    rhs.resize(nb_rows, 0.0);
    // g_array_append_vals(rhs, dbl_values, nb_rows);
    val = lp_interface_get_rhs(RMP, rhs.data());
    lhs_coeff.resize(nb_rows, 0.0);
    // g_array_append_vals(lhs_coeff, dbl_values, nb_rows);
    // id_row = g_array_sized_new(FALSE, FALSE, sizeof(int), nb_rows);
    // g_array_append_vals(id_row, int_values, nb_rows);
    id_row.reserve(nb_rows);
    coeff_row.reserve(nb_rows);

CLEAN:

    if (val) {
        lp_interface_free(&(RMP));
        // g_array_free(pi, TRUE);
        // g_array_free(rhs, TRUE);
        // g_array_free(lhs_coeff, TRUE);
    }
    CC_IFFREE(dbl_values, double);
    CC_IFFREE(int_values, int);
    return val;
}

int NodeData::get_solution_lp_lowerbound() {
    int  val = 0;
    Job* tmp_j = nullptr;

    val = lp_interface_get_nb_cols(RMP, &nb_cols);
    lambda = CC_SAFE_REALLOC(lambda, nb_cols, double);
    std::span lambda_span{lambda, static_cast<size_t>(nb_cols)};
    CCcheck_NULL_2(lambda, "Failed to allocate memory");
    assert(nb_cols == localColPool->len);
    lp_interface_x(RMP, lambda, 0);

    for (int i = 0; i < nb_cols; ++i) {
        ScheduleSet* tmp = ((ScheduleSet*)g_ptr_array_index(localColPool, i));
        if (lambda_span[i]) {
            fmt::print("{}: ", lambda_span[i]);
            g_ptr_array_foreach(tmp->job_list, g_print_machine, NULL);
            fmt::print("\n");
            for (unsigned j = 0; j < tmp->job_list->len; ++j) {
                tmp_j = (Job*)g_ptr_array_index(tmp->job_list, j);
                // printf("%d (%d) ", tmp->num[tmp_j->job],
                //    tmp_j->processing_time);
            }
            fmt::print("\n");
        }
    }

CLEAN:
    return val;
}
