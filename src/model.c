#include <wct.h>

int grab_integer_solution(NodeData* pd, double* x, double tolerance) {
    int          val = 0;
    double       test_incumbent = .0;
    double       incumbent;
    int          i;
    int          tot_weighted = 0;
    int          nb_cols;
    ScheduleSet* tmp_schedule;
    Job*         tmp_j;

    val = wctlp_objval(pd->RMP, &incumbent);
    CCcheck_val_2(val, "wctlp_objval failed");
    val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
    CCcheck_val_2(val, "Failed get nb_cols");

    schedulesets_free(&(pd->bestcolors), &(pd->nb_best));
    pd->bestcolors = CC_SAFE_MALLOC(pd->nb_machines, ScheduleSet);
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    pd->nb_best = 0;

    assert(nb_cols == pd->localColPool->len);
    for (i = 0; i < pd->localColPool->len; ++i) {
        tmp_schedule = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        test_incumbent += x[i];

        if (x[i] >= 1.0 - tolerance) {
            int j = pd->nb_best;
            int k;
            scheduleset_init(pd->bestcolors + j);

            g_ptr_array_set_size(pd->bestcolors[j].job_list,
                                 tmp_schedule->job_list->len);
            for (k = 0; k < tmp_schedule->job_list->len; ++k) {
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

int add_scheduleset_to_rmp(ScheduleSet* set, NodeData* pd) {
    int        val = 0;
    int        row_ind;
    int        var_ind;
    double     cval;
    int        nb_jobs = pd->nb_jobs;
    GPtrArray* members = set->job_list;
    wctlp*     lp = pd->RMP;
    Job*       job;

    val = wctlp_get_nb_cols(lp, &(set->id));
    CCcheck_val_2(val, "Failed to get the number of cols");
    var_ind = set->id;
    val = wctlp_addcol(lp, 0, NULL, NULL,
                       (double)set->total_weighted_completion_time, 0.0,
                       GRB_INFINITY, wctlp_CONT, NULL);
    CCcheck_val_2(val, "Failed to add column to lp")

    for (unsigned i = 0; i < members->len; ++i) {
        job = (Job*)g_ptr_array_index(members, i);
        row_ind = job->job;
        val = wctlp_getcoeff(lp, &row_ind, &var_ind, &cval);
        CCcheck_val_2(val, "Failed wctlp_getcoeff");
        cval += 1.0;
        set->num[job->job] += 1;
        val = wctlp_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);
        CCcheck_val_2(val, "Failed wctlp_chgcoeff");
    }

    row_ind = nb_jobs;
    cval = -1.0;
    val = wctlp_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);
    CCcheck_val_2(val, "Failed wctlp_chgcoeff");

CLEAN:
    return val;
}

void g_add_col_to_lp(gpointer data, gpointer user_data) {
    ScheduleSet* tmp = (ScheduleSet*)data;
    NodeData*    pd = (NodeData*)user_data;
    add_scheduleset_to_rmp(tmp, pd);
}

int build_rmp(NodeData* pd, int construct) {
    int      val = 0;
    Problem* problem = pd->problem;
    int      nb_jobs = problem->nb_jobs;
    int      nb_machines = problem->nb_machines;
    Parms*   parms = &(problem->parms);
    int      nb_row;
    int*     covered = CC_SAFE_MALLOC(nb_jobs, int);
    CCcheck_NULL_2(covered, "Failed to allocate memory to covered");
    fill_int(covered, nb_jobs, 0);
    val = wctlp_init(&(pd->RMP), NULL);
    CCcheck_val_2(val, "wctlp_init failed");

    /**
     * add assignment constraints
     */
    for (int i = 0; i < nb_jobs; i++) {
        val = wctlp_addrow(pd->RMP, 0, (int*)NULL, (double*)NULL,
                           wctlp_EQUAL, 1.0, (char*)NULL);
        CCcheck_val_2(val, "Failed wctlp_addrow");
    }

    /**
     * add number of machines constraint (convexification)
     */
    val = wctlp_addrow(pd->RMP, 0, (int*)NULL, (double*)NULL,
                       wctlp_EQUAL, -(double)nb_machines, (char*)NULL);
    CCcheck_val_2(val, "Failed to add convexification constraint");

    wctlp_get_nb_rows(pd->RMP, &nb_row);

    /** add columns from localColPool */
    g_ptr_array_foreach(pd->localColPool, g_add_col_to_lp, pd);

    /**
     * Some aux variables for column generation
     */
    pd->pi = (double*)CCutil_reallocrus(pd->pi, nb_row * sizeof(double));
    CCcheck_NULL_2(pd->pi, "Failed to allocate memory");
    pd->pi_in = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->pi_in, "Failed to allocate memory");
    fill_dbl(pd->pi_in, nb_row, 0.0);
    pd->eta_in = 0.0;
    pd->pi_out = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->pi_out, "Failed to allocate memory");
    pd->eta_out = 0.0;
    fill_dbl(pd->pi_out, nb_row, 0.0);
    pd->pi_sep = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->pi_sep, "Failed to allocate memory");
    fill_dbl(pd->pi_sep, nb_row, 0.0);
    pd->subgradient_in = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->subgradient_in, "Failed to allocate memory");
    pd->subgradient = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->subgradient, "Failed to allocate memory");
    if (parms->pricing_solver < dp_solver) {
        pd->x_e = CC_SAFE_MALLOC(2 * get_nb_vertices(pd->solver), double);
        CCcheck_NULL_2(pd->x_e, "Failed to allocate memory");
    }
    pd->rhs = CC_SAFE_MALLOC(nb_row, double);
    CCcheck_NULL_2(pd->rhs, "Failed to allocate memory");
    val = wctlp_get_rhs(pd->RMP, pd->rhs);
    CCcheck_val_2(val, "Failed to get RHS");
CLEAN:

    if (val) {
        wctlp_free(&(pd->RMP));
        CC_IFFREE(pd->coeff, double);
        CC_IFFREE(pd->pi, double);
        CC_IFFREE(pd->pi_in, double)
        CC_IFFREE(pd->pi_out, double)
        CC_IFFREE(pd->pi_sep, double)
        CC_IFFREE(pd->subgradient, double)
        CC_IFFREE(pd->subgradient_in, double)
        CC_IFFREE(pd->rhs, double)
    }

    CC_IFFREE(covered, int);
    return val;
}

int get_solution_lp_lowerbound(NodeData* pd) {
    int  val = 0;
    int  nb_cols;
    Job* tmp_j;

    val = wctlp_get_nb_cols(pd->RMP, &nb_cols);
    CCcheck_val_2(val, "Failed in wctlp_get_nb_cols");
    pd->lambda = CC_SAFE_REALLOC(pd->lambda, nb_cols, double);
    CCcheck_NULL_2(pd->lambda, "Failed to allocate memory");
    assert(nb_cols == pd->localColPool->len);
    wctlp_x(pd->RMP, pd->lambda, 0);

    for (unsigned i = 0; i < nb_cols; ++i) {
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
