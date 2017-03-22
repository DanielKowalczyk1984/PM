#include <wct.h>


int print_to_csv(wctproblem *problem) {
    int       val = 0;
    wctdata * pd = &(problem->root_pd);
    wctparms *parms = &(problem->parms);
    FILE *    file = (FILE *)NULL;
    char      filenm[128];
    int       size;
    GDate     date;
    g_date_set_time_t(&date, time(NULL));
    problem->real_time = getRealTime() - problem->real_time;
    CCutil_stop_timer(&(problem->tot_cputime), 0);

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
        case cbfs_conflict_strategy:
            sprintf(filenm, "WCT_CONFLICT_%d_%d.csv", pd->nmachines, pd->njobs);
            break;

        case ahv_strategy:
        case cbfs_ahv_strategy:
            sprintf(filenm, "WCT_AHV_%d_%d.csv", pd->nmachines, pd->njobs);
            break;
    }

    file = fopen(filenm, "a+");

    if (file == NULL) {
        printf("We couldn't open %s in %s at line %d\n", filenm, __FILE__,
               __LINE__);
        val = 1;
        goto CLEAN;
    }

    fseek(file, 0, SEEK_END);
    size = ftell(file);

    if (size == 0) {
        fprintf(file,
                "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%"
                "s,%s\n",
                "NameInstance", "tot_real_time", "tot_cputime", "tot_lb",
                "tot_lb_root", "tot_lb_lp", "tot_branch_and_bound",
                "tot_scatter_search", "tot_build_dd", "tot_pricing",
                "rel_error", "status", "global_lower_bound",
                "global_upper_bound", "first_lower_bound", "first_upper_bound",
                "first_rel_error", "solved_at_root", "nb_explored_nodes",
                "nb_generated_col", "nb_generated_col_root", "date");
    }

    fprintf(file,
            "%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%f,%d,%d,%d,%d,%u/"
            "%u/%u\n",
            pd->pname, problem->real_time, problem->tot_cputime.cum_zeit,
            problem->tot_lb.cum_zeit, problem->tot_lb_lp_root.cum_zeit,
            problem->tot_lb_lp.cum_zeit, problem->tot_branch_and_bound.cum_zeit,
            problem->tot_scatter_search.cum_zeit,
            problem->tot_build_dd.cum_zeit, problem->tot_pricing.cum_zeit,
            problem->rel_error, problem->status, problem->global_lower_bound,
            problem->global_upper_bound, problem->first_lower_bound,
            problem->first_upper_bound, problem->first_rel_error,
            problem->found, problem->nb_explored_nodes,
            problem->nb_generated_col, problem->nb_generated_col_root, date.day,
            date.month, date.year);
    fclose(file);
CLEAN:
    return val;
}

int print_to_screen(wctproblem *problem) {
    int val = 0;

    switch (problem->status) {
        case no_sol:
            printf(
                "We didn't decide if this instance is feasible or "
                "infeasible\n");
            break;

        case feasible:
        case lp_feasible:
        case meta_heur:
            printf("A suboptimal schedule with relative error %f is found.\n",
                   (double)(problem->global_upper_bound -
                            problem->global_lower_bound) /
                       (problem->global_lower_bound));
            break;

        case optimal:
            printf("The optimal schedule is found.\n");
            break;
    }

    printf(
        "Compute_schedule took %f seconds(tot_scatter_search %f, "
        "tot_branch_and_bound %f, tot_lb_lp_root %f, tot_lb_lp %f, tot_lb %f, "
        "tot_pricing %f, tot_build_dd %f) and %f seconds in real time\n",
        problem->tot_cputime.cum_zeit, problem->tot_scatter_search.cum_zeit,
        problem->tot_branch_and_bound.cum_zeit,
        problem->tot_lb_lp_root.cum_zeit, problem->tot_lb_lp.cum_zeit,
        problem->tot_lb.cum_zeit, problem->tot_pricing.cum_zeit,
        problem->tot_build_dd.cum_zeit, problem->real_time);
    return val;
}