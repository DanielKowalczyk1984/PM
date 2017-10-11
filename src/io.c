#include <wct.h>

static int get_problem_name(char *pname, const char *efname) {
    int         rval = 0;
    int         len = 0;
    const char *fname = strrchr(efname, '/');
    const char *lastdot = strrchr(efname, '.');

    if (!fname) {
        fname = efname;
    } else {
        fname++;
    }

    if (lastdot) {
        len = lastdot - fname + 1;
    } else {
        len = strlen(fname);
    }

    if (snprintf(pname, len, "%s", fname) < 0) {
        rval = 1;
    }

    printf("Extracted problem name %s\n", pname);
    return rval;
}

int read_problem(wctproblem *problem) {
    int         val = 0;
    int         nbjobs = 0;
    int         curduration, curduedate, curweight, curjob;
    Job *       _jobarray = (Job *)NULL;
    Job * tmp_j;
    char        buf[256], *p;
    int         bufsize;
    const char *delim = " \n\t";
    char *      data = (char *)NULL;
    char *      buf2 = (char *)NULL;
    wctdata *   pd;
    wctparms *  parms;
    parms = &(problem->parms);
    pd = &(problem->root_pd);
    FILE *in = fopen(parms->jobfile, "r");
    curjob = 0;

    if (in != (FILE *)NULL) {
        get_problem_name(pd->pname, parms->jobfile);

        if (fgets(buf, 254, in) != NULL) {
            p = buf;
            data = strtok(p, delim);
            sscanf(data, "%d", &nbjobs);
            bufsize = 3 * nbjobs * (2 + (int)ceil(log((double)nbjobs + 10)));
            buf2 = (char *)CC_SAFE_MALLOC(bufsize, char);
            CCcheck_NULL_2(buf2, "Failed to allocate buf2");
        } else {
            val = 1;
            goto CLEAN;
        }

        while (fgets(buf2, bufsize, in) != (char *)NULL) {
            p = buf2;
            data = strtok(p, delim);
            sscanf(data, "%d", &curduration);
            data = strtok(NULL, delim);
            sscanf(data, "%d", &curduedate);
            data = strtok(NULL, delim);
            sscanf(data, "%d", &curweight);
            curduedate = curduedate/parms->nmachines;
            tmp_j = job_alloc(&curduration, &curweight, &curduedate);
            g_ptr_array_add(problem->g_job_array, tmp_j);

            if (tmp_j->processingime > tmp_j->duetime) {
                problem->off += tmp_j->weight * (tmp_j->processingime - tmp_j->duetime);
                tmp_j->duetime = tmp_j->processingime;
            }

            tmp_j->job = curjob;
            curjob++;
        }

        problem->njobs = pd->njobs = nbjobs;
        problem->nmachines = parms->nmachines;
    } else {
        fprintf(stderr, "Unable to open file %s\n", parms->jobfile);
        val = 1;
        goto CLEAN;
    }

CLEAN:

    if (val) {
        CC_IFFREE(_jobarray, Job);
    }

    CC_IFFREE(buf2, char);

    if (in) {
        fclose(in);
    }

    return val;
}

int print_to_csv(wctproblem *problem) {
    int       val = 0;
    wctdata * pd = &(problem->root_pd);
    wctparms *parms = &(problem->parms);
    FILE *    file = (FILE *)NULL;
    char      filenm[128];
    int       size;
    GDate     date;
    g_date_set_time_t(&date, time(NULL));
    problem->real_time_total = getRealTime() - problem->real_time_total;
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
            pd->pname, problem->real_time_total, problem->tot_cputime.cum_zeit,
            problem->tot_lb.cum_zeit, problem->tot_lb_root.cum_zeit,
            problem->tot_lb.cum_zeit, problem->tot_branch_and_bound.cum_zeit,
            problem->tot_heuristic.cum_zeit,
            problem->tot_build_dd.cum_zeit, problem->tot_pricing.cum_zeit,
            problem->rel_error, problem->status, problem->global_lower_bound,
            problem->global_upper_bound, problem->root_lower_bound,
            problem->root_upper_bound, problem->root_rel_error,
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
        problem->tot_cputime.cum_zeit, problem->tot_heuristic.cum_zeit,
        problem->tot_branch_and_bound.cum_zeit,
        problem->tot_lb_root.cum_zeit, problem->tot_lb.cum_zeit,
        problem->tot_lb.cum_zeit, problem->tot_pricing.cum_zeit,
        problem->tot_build_dd.cum_zeit, problem->real_time_total);
    return val;
}

/** Printing sizes of ZDD */
int print_size_to_csv(wctproblem *problem, wctdata *pd) {
    int       val = 0;
    int       size;
    wctdata * root_node = &(problem->root_pd);
    wctparms *parms = &(problem->parms);
    char      filenm[128];
    FILE *    file = (FILE *)NULL;
    GDate     date;
    g_date_set_time_t(&date, time(NULL));

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
        case cbfs_conflict_strategy:
            sprintf(filenm, "CONFLICT_%d_%d.csv", pd->nmachines, pd->njobs);
            break;

        case ahv_strategy:
        case cbfs_ahv_strategy:
            sprintf(filenm, "AHV_%d_%d.csv", pd->nmachines, pd->njobs);
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

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
            if (size == 0) {
                fprintf(file, "%s;%s;%s;%s;%s;%s;%s;%s\n", "NameInstance",
                        "depth", "size", "nb_same", "nb_diff", "v1", "v2",
                        "date");
            }

            fprintf(file, "%s;%d;%zu;%d;%d;%d;%d;%u/%u/%u\n", root_node->pname,
                    pd->depth, get_datasize(pd->solver), pd->ecount_same,
                    pd->ecount_differ, pd->v1->job, pd->v2->job, date.day, date.month,
                    date.year);
            break;

        case ahv_strategy:
            if (size == 0) {
                fprintf(file, "%s;%s;%s;%s;%s;%s\n", "NameInstance", "depth",
                        "size", "branch_job", "completion_time", "data");
            }

            fprintf(file, "%s;%d;%zu;%d;%d;%u/%u/%u\n", root_node->pname,
                    pd->depth, get_datasize(pd->solver), pd->branch_job,
                    pd->completiontime, date.day, date.month, date.year);
            break;
    }

    fclose(file);
CLEAN:
    return val;
}
