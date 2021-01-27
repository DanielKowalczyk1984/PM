#include <wct.h>
#include "gurobi_c.h"
#include "pricingstabilizationwrapper.h"
#include "scheduleset.h"
#include "solver.h"
#include "util.h"
#include "wctprivate.h"

int debug = 0;

/*Information about debug*/
int dbg_lvl() {
    return debug;
}
void set_dbg_lvl(int dbglvl) {
    debug = dbglvl;
}

/*Functions for initialization of the problem and freeing the problem*/
void problem_init(Problem* problem) {
    problem->root_pd = CC_SAFE_MALLOC(1, NodeData);
    nodedata_init(problem->root_pd, problem);
    parms_init(&(problem->parms));
    statistics_init(&(problem->stat));
    problem->tree = (BranchBoundTree*)NULL;
    /** Job data */
    problem->g_job_array = g_ptr_array_new_with_free_func(g_job_free);
    problem->intervals = g_ptr_array_new_with_free_func(g_interval_free);
    problem->opt_sol = (Solution*)NULL;
    /** Job summary */
    problem->nb_jobs = 0;
    problem->p_sum = 0;
    problem->pmax = 0;
    problem->pmin = INT_MAX;
    problem->dmax = INT_MIN;
    problem->dmin = INT_MAX;
    problem->off = 0;
    problem->H_min = 0;
    problem->H_max = INT_MAX;
    /*B&B info*/
    problem->nb_data_nodes = 0;
    problem->global_upper_bound = INT_MAX;
    problem->global_lower_bound = 0;
    problem->rel_error = DBL_MAX;
    problem->root_lower_bound = 0;
    problem->root_upper_bound = INT_MAX;
    problem->root_rel_error = DBL_MAX;
    problem->status = no_sol;
    problem->br_heap_a = (BinomialHeap*)NULL;
    /*data of the problem*/
    set_id_and_name(problem->root_pd, 0, "root_node");
    problem->nb_data_nodes++;
    /*parms of the problem*/
    /** statistics of the problem */
    /*heap initialization*/
    problem->unexplored_states = g_ptr_array_new();
    problem->non_empty_level_pqs = g_queue_new();
    problem->last_explored = -1;
    problem->found = 0;
    /** initialize colPool */
    problem->ColPool = g_ptr_array_new_with_free_func(g_scheduleset_free);
}

void problem_free(Problem* problem) {
    /*free the parameters*/
    parms_free(&(problem->parms));
    delete_branch_bound_tree(problem->tree);
    // nodedata_free(problem->root_pd);

    /*free the heap*/
    if (problem->br_heap_a != (BinomialHeap*)NULL) {
        binomial_heap_free(problem->br_heap_a);
    }

    for (unsigned int i = 0; i < problem->unexplored_states->len; ++i) {
        if ((problem->unexplored_states->pdata[i]) != (BinomialHeap*)NULL) {
            binomial_heap_free(
                (BinomialHeap*)problem->unexplored_states->pdata[i]);
        }
    }

    g_ptr_array_free(problem->g_job_array, TRUE);
    g_ptr_array_free(problem->unexplored_states, TRUE);
    g_ptr_array_free(problem->ColPool, TRUE);
    g_ptr_array_free(problem->intervals, TRUE);
    g_queue_free(problem->non_empty_level_pqs);
    solution_free(&(problem->opt_sol));
}

/*Functions for initialization and free the data*/
void nodedata_init(NodeData* pd, Problem* prob) {
    /*Initialization B&B data*/
    pd->id = -1;
    pd->depth = 0;
    pd->status = initialized;
    sprintf(pd->pname, "temporary");
    /*Initialization node instance data*/
    pd->nb_jobs = 0;
    pd->H_max = 0;
    pd->H_min = 0;
    pd->off = prob->off;
    pd->ordered_jobs = g_ptr_array_new_with_free_func(g_free);
    pd->jobarray = (GPtrArray*)NULL;
    /** Initialization data */
    pd->upper_bound = INT_MAX;
    pd->lower_bound = 0;
    pd->LP_lower_bound = 0.0;
    pd->LP_lower_min = DBL_MAX;
    pd->rhs = (GArray*)NULL;
    /*Initialization  of the LP*/
    pd->RMP = (wctlp*)NULL;
    pd->lambda = (double*)NULL;
    // pd->x_e = (double*)NULL;
    pd->pi = (GArray*)NULL;
    pd->slack = (GArray*)NULL;
    pd->lhs_coeff = (GArray*)NULL;
    pd->id_row = (GArray*)NULL;
    pd->coeff_row = (GArray*)NULL;
    pd->nb_rows = 0;
    pd->nb_cols = 0;
    // init info cut generation
    pd->max_nb_cuts = NB_CUTS;
    pd->id_convex_constraint = 0;
    pd->id_assignment_constraint = 0;
    pd->id_valid_cuts = 0;
    pd->id_art_var_assignment = 0;
    pd->id_art_var_convex = 0;
    pd->id_art_var_cuts = 0;
    pd->id_pseudo_schedules = 0;

    /**init stab data */
    pd->update = 1;
    pd->iterations = 0;
    /*Initialization pricing_problem*/
    pd->solver = (PricerSolver*)NULL;
    pd->nb_non_improvements = 0;
    pd->zero_count = 0;
    pd->best_schedule = g_ptr_array_new_with_free_func(g_scheduleset_free);
    /**Column schedules */
    pd->localColPool = g_ptr_array_new_with_free_func(g_scheduleset_free);
    pd->column_status = (int*)NULL;
    /*Initialization max and retirement age*/
    pd->maxiterations = NB_CG_ITERATIONS;
    pd->retirementage = 100;
    /*initialization of branches*/
    pd->branch_job = -1;
    pd->less = -1;
    pd->completiontime = 0;
    pd->parent = (NodeData*)NULL;
    pd->parms = &(prob->parms);
}

void nodedata_init_null(NodeData* pd) {
    /*Initialization B&B data*/
    pd->id = -1;
    pd->depth = 0;
    pd->status = initialized;
    sprintf(pd->pname, "temporary");
    /*Initialization node instance data*/
    pd->nb_jobs = 0;
    pd->H_max = 0;
    pd->H_min = 0;
    pd->off = 0;
    pd->ordered_jobs = (GPtrArray*)NULL;
    pd->jobarray = (GPtrArray*)NULL;
    /** Initialization data */
    pd->upper_bound = INT_MAX;
    pd->lower_bound = 0;
    pd->LP_lower_bound = 0.0;
    pd->LP_lower_min = DBL_MAX;
    pd->rhs = (GArray*)NULL;
    /*Initialization  of the LP*/
    pd->RMP = (wctlp*)NULL;
    pd->lambda = (double*)NULL;
    // pd->x_e = (double*)NULL;
    pd->pi = (GArray*)NULL;
    pd->slack = (GArray*)NULL;
    pd->lhs_coeff = (GArray*)NULL;
    pd->id_row = (GArray*)NULL;
    pd->coeff_row = (GArray*)NULL;
    pd->nb_rows = 0;
    pd->nb_cols = 0;
    // init info cut generation
    pd->max_nb_cuts = NB_CUTS;
    pd->id_convex_constraint = 0;
    pd->id_assignment_constraint = 0;
    pd->id_valid_cuts = 0;
    pd->id_art_var_assignment = 0;
    pd->id_art_var_convex = 0;
    pd->id_art_var_cuts = 0;
    pd->id_pseudo_schedules = 0;

    /**init stab data */
    pd->update = 1;
    pd->iterations = 0;
    /*Initialization pricing_problem*/
    pd->solver = (PricerSolver*)NULL;
    pd->nb_non_improvements = 0;
    pd->zero_count = 0;
    // pd->bestcolors = (ScheduleSet*)NULL;
    pd->best_schedule = (GPtrArray*)NULL;
    // pd->nb_best = 0;
    /**Column schedules */
    pd->localColPool = (GPtrArray*)NULL;
    pd->column_status = (int*)NULL;
    /*Initialization max and retirement age*/
    pd->maxiterations = NB_CG_ITERATIONS;
    pd->retirementage = 100;
    /*initialization of branches*/
    pd->branch_job = -1;
    pd->parent = (NodeData*)NULL;
    pd->completiontime = 0;
    pd->less = -1;
    pd->parms = (Parms*)NULL;
}

NodeData* new_node_data(NodeData* pd) {
    NodeData* aux = CC_SAFE_MALLOC(1, NodeData);
    nodedata_init_null(aux);

    aux->parms = pd->parms;
    aux->stat = pd->stat;
    aux->opt_sol = pd->opt_sol;
    aux->depth = pd->depth + 1;

    /** Instance copy */
    aux->jobarray = pd->jobarray;
    aux->nb_jobs = pd->nb_jobs;
    aux->nb_machines = pd->nb_machines;
    aux->H_max = pd->H_max;
    aux->H_min = pd->H_min;
    aux->off = pd->off;

    /** copy info about intervals */
    aux->ordered_jobs =
        g_ptr_array_copy(pd->ordered_jobs, g_copy_interval_pair, NULL);

    /** info about RMP */

    /** Ids for model */
    aux->max_nb_cuts = pd->max_nb_cuts;
    aux->id_convex_constraint = pd->id_convex_constraint;
    aux->id_assignment_constraint = pd->id_assignment_constraint;
    aux->id_valid_cuts = pd->id_valid_cuts;

    aux->id_art_var_convex = pd->id_art_var_convex;
    aux->id_art_var_assignment = pd->id_art_var_assignment;
    aux->id_art_var_cuts = pd->id_art_var_cuts;
    aux->id_next_var_cuts = pd->id_next_var_cuts;
    aux->id_pseudo_schedules = pd->id_pseudo_schedules;

    /** copy info about solver */

    aux->solver = copy_pricer_solver(pd->solver, aux->ordered_jobs, aux->parms);

    aux->localColPool =
        g_ptr_array_copy(pd->localColPool, g_copy_scheduleset, &(pd->nb_jobs));

    aux->lower_bound = pd->lower_bound;
    aux->upper_bound = pd->upper_bound;

    aux->LP_lower_bound = pd->LP_lower_bound;
    aux->LP_lower_bound_dual = pd->LP_lower_bound_dual;
    aux->LP_lower_bound_BB = pd->LP_lower_bound_BB;
    aux->LP_lower_min = pd->LP_lower_min;

    aux->iterations = 0;
    aux->nb_non_improvements = 0;

    aux->solver_stab = new_pricing_stabilization(aux->solver, aux->parms);
    /** copy info about best_schedule */
    aux->best_schedule =
        g_ptr_array_copy(pd->best_schedule, g_copy_scheduleset, NULL);

    aux->maxiterations = pd->maxiterations;
    aux->retirementage = pd->retirementage;

    return aux;
}

void lp_node_data_free(NodeData* pd) {
    /**
     * free all the gurobi data associated with the LP relaxation
     */
    if (pd->RMP) {
        lp_interface_free(&(pd->RMP));
    }

    /**
     * free all the data associated with the LP
     */
    if (pd->pi) {
        g_array_free(pd->pi, TRUE);
    }
    if (pd->slack) {
        g_array_free(pd->slack, TRUE);
    }
    CC_IFFREE(pd->lambda, double);
    if (pd->rhs) {
        g_array_free(pd->rhs, TRUE);
    }

    if (pd->lhs_coeff) {
        g_array_free(pd->lhs_coeff, TRUE);
    }

    if (pd->id_row) {
        g_array_free(pd->id_row, TRUE);
    }
    if (pd->coeff_row) {
        g_array_free(pd->coeff_row, TRUE);
    }
    CC_IFFREE(pd->column_status, int);

    /**
     * free all the schedules from the localColPool
     */
    g_ptr_array_free(pd->localColPool, TRUE);
    pd->nb_rows = 0;
    pd->nb_cols = 0;
    pd->max_nb_cuts = NB_CUTS;
    pd->id_convex_constraint = 0;
    pd->id_assignment_constraint = 0;
    pd->id_valid_cuts = 0;
    pd->id_art_var_assignment = 0;
    pd->id_art_var_convex = 0;
    pd->id_art_var_cuts = 0;
    pd->id_pseudo_schedules = 0;
}

void temporary_data_free(NodeData* pd) {
    lp_node_data_free(pd);
    // g_ptr_array_free(pd->localColPool, TRUE);
    if (pd->solver) {
        freeSolver(pd->solver);
        pd->solver = (PricerSolver*)NULL;
    }
    if (pd->solver_stab) {
        delete_pricing_stabilization(pd->solver_stab);
    }
}

void nodedata_free(NodeData* pd) {
    temporary_data_free(pd);

    g_ptr_array_free(pd->ordered_jobs, TRUE);
    g_ptr_array_free(pd->best_schedule, TRUE);
    CC_IFFREE(pd, NodeData);
}

int set_id_and_name(NodeData* pd, int id, const char* fname) {
    int val = 0;
    int s_val = 0;
    CCcheck_NULL_2(pd, "np memory was allocated to pd");
    pd->id = id;
    s_val = snprintf(pd->pname, MAX_PNAME_LEN, "%s", fname);

    if (s_val < 0 || MAX_PNAME_LEN <= s_val) {
        val = 1;
        CCcheck_val(val, "Failed to write pname")
    }

CLEAN:
    return val;
}
static void scheduleset_unify(GPtrArray* array) {
    int          i = 0;
    int          it = 1;
    int          first_del = -1;
    int          last_del = -1;
    int          nb_col = array->len;
    ScheduleSet* temp = (ScheduleSet*)NULL;
    ScheduleSet* prev = (ScheduleSet*)NULL;
    g_ptr_array_sort(array, g_scheduleset_less);

    if (!(array->len)) {
        return;
    }

    prev = (ScheduleSet*)g_ptr_array_index(array, 0);
    /* Find first non-empty set */
    for (i = 1; i < nb_col; ++i) {
        temp = (ScheduleSet*)g_ptr_array_index(array, it);
        if (scheduleset_less(prev, temp)) {
            if (first_del != -1) {
                /** Delete recently found deletion range.*/
                g_ptr_array_remove_range(array, first_del,
                                         last_del - first_del + 1);
                it = it - (last_del - first_del);
                first_del = last_del = -1;
            } else {
                it++;
            }
            prev = temp;
        } else {
            if (first_del == -1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            prev = temp;
            it++;
        }
    }

    if (first_del != -1) {
        g_ptr_array_remove_range(array, first_del, last_del - first_del + 1);
    }
}

int prune_duplicated_sets(NodeData* pd) {
    int val = 0;
    scheduleset_unify(pd->localColPool);

    if (dbg_lvl() > 1) {
        for (guint i = 0; i < pd->localColPool->len; ++i) {
            ScheduleSet* tmp =
                (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
            GPtrArray* tmp_a = tmp->job_list;

            printf("TRANSSORT SET ");

            for (guint j = 0; j < tmp_a->len; ++j) {
                Job* tmp_j = (Job*)g_ptr_array_index(tmp_a, j);
                printf(" %d", tmp_j->job);
            }

            printf("\n");
        }
    }

    return val;
}

// int compute_schedule(Problem* problem) {
//     int         val = 0;
//     NodeData*   root_pd = problem->root_pd;
//     Parms*      parms = &(problem->parms);
//     Statistics* statistics = &(problem->stat);
//     problem->mult_key = 1.0;
//     problem->root_upper_bound = problem->global_upper_bound;
//     problem->root_lower_bound = problem->global_lower_bound;
//     problem->root_rel_error =
//         (double)(problem->global_upper_bound - problem->global_lower_bound) /
//         ((double)problem->global_lower_bound + 0.00001);
//     prune_duplicated_sets(root_pd);
//     init_BB_tree(problem);

//     if (root_pd->status >= LP_bound_computed) {
//         val = prefill_heap(root_pd, problem);
//         CCcheck_val(val, "Failed in prefill_heap");
//     } else {
//         CCutil_start_timer(&(statistics->tot_lb_root));
//         val = compute_lower_bound(root_pd);
//         CCcheck_val_2(val, "Failed in compute_lower_bound");

//         if (root_pd->lower_bound > problem->global_lower_bound) {
//             problem->global_lower_bound = root_pd->lower_bound;
//             problem->root_lower_bound = root_pd->lower_bound;
//             problem->root_rel_error =
//                 (double)(problem->root_upper_bound -
//                          problem->root_lower_bound) /
//                 ((double)problem->root_lower_bound + 0.00001);
//         }

//         CCcheck_val_2(val, "Failed in compute_lower_bound");
//         statistics->nb_generated_col_root = statistics->nb_generated_col;
//         CCutil_stop_timer(&(statistics->tot_lb_root), 0);

//         switch (parms->bb_branch_strategy) {
//             case conflict_strategy:
//             case ahv_strategy:
//                 val = insert_into_branching_heap(root_pd, problem);
//                 CCcheck_val_2(val, "insert_into_branching_heap failed");
//                 break;

//             case cbfs_conflict_strategy:
//             case cbfs_ahv_strategy:
//                 insert_node_for_exploration(root_pd, problem);
//                 break;
//         }
//     }

//     printf("GUB = %d, GLB = %d\n", problem->global_upper_bound,
//            problem->global_lower_bound);

//     if (problem->global_lower_bound != problem->global_upper_bound) {
//         CCutil_start_resume_time(&(statistics->tot_branch_and_bound));

//         switch (parms->bb_branch_strategy) {
//             case conflict_strategy:
//                 val = sequential_branching_conflict(problem);
//                 CCcheck_val(val, "Failed in sequential_branching_conflict");
//                 break;

//             case cbfs_conflict_strategy:
//                 val = sequential_cbfs_branch_and_bound_conflict(problem);
//                 CCcheck_val_2(val, "Failed in CBFS conflict branching");
//                 break;
//         }

//         CCutil_stop_timer(&(statistics->tot_branch_and_bound), 0);
//         printf("Compute schedule finished with LB %d and UB %d\n",
//                root_pd->lower_bound, problem->global_upper_bound);
//     } else {
//         problem->found = 1;
//     }

//     if (root_pd->lower_bound == problem->global_upper_bound) {
//         problem->global_lower_bound = root_pd->lower_bound;
//         problem->rel_error = (double)(problem->global_upper_bound -
//                                       problem->global_lower_bound) /
//                              ((double)problem->global_lower_bound);
//         problem->status = optimal;
//         printf("The optimal schedule is given by:\n");
//         // print_schedule(root_pd->bestcolors, root_pd->nb_best);
//         printf("with total weighted completion time %d\n",
//                root_pd->upper_bound);
//     } else {
//         problem->global_lower_bound = root_pd->lower_bound;
//         problem->rel_error = (double)(problem->global_upper_bound -
//                                       problem->global_lower_bound) /
//                              ((double)problem->global_lower_bound);
//         problem->status = meta_heuristic;
//         problem->global_lower_bound = root_pd->lower_bound;
//         printf("The suboptimal schedule is given by:\n");
//         // print_schedule(root_pd->bestcolors, root_pd->nb_best);
//         printf("with total weighted completion time\n");
//     }

//     children_data_free(problem->root_pd);
// CLEAN:
//     return val;
// }

int add_solution_to_colpool(Solution* sol, NodeData* pd) {
    int val = 0;

    for (int i = 0; i < sol->nb_machines; ++i) {
        GPtrArray*   machine = sol->part[i].machine;
        ScheduleSet* tmp = scheduleset_from_solution(machine, pd->nb_jobs);
        CCcheck_NULL_2(tmp, "Failed to allocate memory");
        tmp->id = pd->localColPool->len;
        g_ptr_array_add(pd->localColPool, tmp);
    }

CLEAN:
    return val;
}

int add_solution_to_colpool_and_lp(Solution* sol, NodeData* pd) {
    int          val = 0;
    ScheduleSet* tmp = (ScheduleSet*)NULL;

    for (int i = 0; i < sol->nb_machines; ++i) {
        GPtrArray* machine = sol->part[i].machine;
        tmp = scheduleset_from_solution(machine, pd->nb_jobs);
        CCcheck_NULL_2(tmp, "Failed to allocate memory");
        g_ptr_array_add(pd->localColPool, tmp);
    }

    for (unsigned i = 0; i < pd->localColPool->len; ++i) {
        tmp = (ScheduleSet*)g_ptr_array_index(pd->localColPool, i);
        add_scheduleset_to_rmp(tmp, pd);
    }

CLEAN:
    return val;
}
