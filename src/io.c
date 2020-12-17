#include <unistd.h>
#include <wct.h>
static int get_problem_name(char* pname, const char* end_file_name) {
    int         rval = 0;
    int         len = 0;
    const char* fname = strrchr(end_file_name, '/');
    const char* lastdot = strrchr(end_file_name, '.');

    if (!fname) {
        fname = end_file_name;
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

int read_problem(Problem* problem) {
    int         val = 0;
    int         nb_jobs = 0;
    int         curduration, curduedate, curweight, curjob;
    Job*        _jobarray = (Job*)NULL;
    Job*        tmp_j;
    char        buf[256], *p;
    int         bufsize;
    const char* delim = " \n\t";
    char*       data = (char*)NULL;
    char*       buf2 = (char*)NULL;
    NodeData*   pd;
    Parms*      parms;
    parms = &(problem->parms);
    pd = (problem->root_pd);
    FILE* in = fopen(parms->jobfile, "r");
    curjob = 0;

    if (in != (FILE*)NULL) {
        get_problem_name(pd->pname, parms->jobfile);

        if (fgets(buf, 254, in) != NULL) {
            p = buf;
            data = strtok(p, delim);
            sscanf(data, "%d", &nb_jobs);
            bufsize = 3 * nb_jobs * (2 + (int)ceil(log((double)nb_jobs + 10)));
            buf2 = (char*)CC_SAFE_MALLOC(bufsize, char);
            CCcheck_NULL_2(buf2, "Failed to allocate buf2");
        } else {
            val = 1;
            goto CLEAN;
        }

        while (fgets(buf2, bufsize, in) != (char*)NULL) {
            p = buf2;
            sscanf(p, "%d %d %d", &curduration, &curduedate, &curweight);
            curduedate = curduedate / parms->nb_machines;
            tmp_j = job_alloc(&curduration, &curweight, &curduedate);
            g_ptr_array_add(problem->g_job_array, tmp_j);

            if (tmp_j->processing_time > tmp_j->due_time) {
                problem->off +=
                    tmp_j->weight * (tmp_j->processing_time - tmp_j->due_time);
                tmp_j->due_time = tmp_j->processing_time;
            }

            tmp_j->job = curjob;
            curjob++;
        }

        problem->nb_jobs = pd->nb_jobs = parms->nb_jobs = nb_jobs;
        problem->nb_machines = pd->nb_machines = parms->nb_machines;
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

int print_to_csv(Problem* problem) {
    int         val = 0;
    NodeData*   pd = problem->root_pd;
    Parms*      parms = &(problem->parms);
    Statistics* statistics = &(problem->stat);
    FILE*       file = (FILE*)NULL;
    char*       file_name = CC_SAFE_MALLOC(128, char);
    GDate       date;
    g_date_set_time_t(&date, time(NULL));
    statistics->real_time_total = getRealTime() - statistics->real_time_total;
    CCutil_stop_timer(&(statistics->tot_cputime), 0);

    sprintf(file_name, "CG_overall_%d%02d%02d.csv", date.year, date.month,
            date.day);

    if (access(file_name, F_OK) != -1) {
        file = fopen(file_name, "a");
    } else {
        file = fopen(file_name, "w");
        if (file == NULL) {
            printf("We couldn't open %s in %s at line %d\n", "CG_overall.csv",
                   __FILE__, __LINE__);
            val = 1;
            goto CLEAN;
        }
        fprintf(
            file,
            "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%"
            "s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
            "NameInstance", "tot_real_time", "tot_cputime", "tot_lb",
            "tot_lb_root", "tot_heuristic", "tot_build_dd", "tot_pricing",
            statistics->tot_reduce_cost_fixing.name, "rel_error",
            "global_lower_bound", "global_upper_bound", "first_rel_error",
            "nb_generated_col", "date", "nb_iterations_rvnd", "stabilization",
            "alpha", "pricing_solver", "n", "m", "first_size_graph",
            "size_after_reduced_cost", "mip_nb_vars", "mip_nb_constr",
            "mip_obj_bound", "mip_obj_bound_lp", "mip_rel_gap", "mip_run_time",
            "mip_status", "mip_nb_iter_simplex", "mip_nb_nodes");
    }

    for (int i = MIP_Attr_Run_Time; i <= MIP_Attr_Nb_Nodes && parms->mip_solver;
         i++) {
        get_mip_statistics(pd, i);
    }

    fprintf(file,
            "%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%d,%u/"
            "%u/%u,%d,%d,%f,%d,%d,%d,%lu,%lu,%d,%d,%f,%f,%f,%f,%d,%f,%f\n",
            pd->pname, statistics->real_time_total,
            statistics->tot_cputime.cum_zeit, statistics->tot_lb.cum_zeit,
            statistics->tot_lb_root.cum_zeit,
            statistics->tot_heuristic.cum_zeit,
            statistics->tot_build_dd.cum_zeit, statistics->tot_pricing.cum_zeit,
            statistics->tot_reduce_cost_fixing.cum_zeit, statistics->rel_error,
            statistics->global_lower_bound, statistics->global_upper_bound,
            statistics->root_rel_error, statistics->nb_generated_col, date.day,
            date.month, date.year, parms->nb_iterations_rvnd,
            parms->stab_technique, parms->alpha, parms->pricing_solver,
            problem->nb_jobs, parms->nb_machines, statistics->first_size_graph,
            statistics->size_graph_after_reduced_cost_fixing,
            statistics->mip_nb_vars, statistics->mip_nb_constr,
            statistics->mip_obj_bound, statistics->mip_obj_bound_lp,
            statistics->mip_rel_gap, statistics->mip_run_time,
            statistics->mip_status, statistics->mip_nb_iter_simplex,
            statistics->mip_nb_nodes);
    fclose(file);
CLEAN:
    CC_FREE(file_name, char);
    return val;
}

int print_to_screen(Problem* problem) {
    int         val = 0;
    Statistics* statistics = &(problem->stat);

    switch (problem->status) {
        case no_sol:
            printf(
                "We didn't decide if this instance is feasible or "
                "infeasible\n");
            break;

        case feasible:
        case lp_feasible:
        case meta_heuristic:
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
        statistics->tot_cputime.cum_zeit, statistics->tot_heuristic.cum_zeit,
        statistics->tot_branch_and_bound.cum_zeit,
        statistics->tot_lb_root.cum_zeit, statistics->tot_lb.cum_zeit,
        statistics->tot_lb.cum_zeit, statistics->tot_pricing.cum_zeit,
        statistics->tot_build_dd.cum_zeit, statistics->real_time_total);
    return val;
}

/** Printing sizes of ZDD */
int print_size_to_csv(Problem* problem, NodeData* pd) {
    int    val = 0;
    int    size;
    Parms* parms = &(problem->parms);
    char   file_name[128];
    FILE*  file = (FILE*)NULL;
    GDate  date;
    g_date_set_time_t(&date, time(NULL));

    sprintf(file_name, "overall_size.csv");

    file = fopen(file_name, "a+");

    if (file == NULL) {
        printf("We couldn't open %s in %s at line %d\n", file_name, __FILE__,
               __LINE__);
        val = 1;
        goto CLEAN;
    }

    fseek(file, 0, SEEK_END);
    size = ftell(file);

    if (size == 0) {
        fprintf(file, "%s,%s,%s,%s,%s\n", "NameInstance", "date", "solver",
                "size", "depth");
    }

    fprintf(file, "%s,%u/%u/%u,%d,%lu,%d\n", pd->pname, date.day, date.month,
            date.year, parms->pricing_solver, get_nb_vertices(pd->solver),
            pd->depth);

    fclose(file);
CLEAN:
    return val;
}
