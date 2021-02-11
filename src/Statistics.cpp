#include "Statistics.h"

void statistics_init(Statistics* stat) {
    stat->global_upper_bound = INT_MAX;
    stat->global_lower_bound = 0;
    stat->rel_error = DBL_MAX;
    stat->root_lower_bound = 0.0;
    stat->root_upper_bound = INT_MAX;
    stat->root_rel_error = DBL_MAX;
    stat->nb_explored_nodes = 0;
    stat->nb_generated_nodes = 0;
    stat->nb_generated_col = 0;
    stat->nb_generated_col_root = 0;

    CCutil_init_timer(&(stat->tot_build_dd), "tot_build_dd");
    CCutil_init_timer(&(stat->tot_cputime), "tot_cputime");
    CCutil_init_timer(&(stat->tot_branch_and_bound), "tot_branch_and_bound");
    CCutil_init_timer(&(stat->tot_strong_branching), "tot_strong_branching");
    CCutil_init_timer(&(stat->tot_lb_root), "tot_lb_root");
    CCutil_init_timer(&(stat->tot_lb), "tot_lb");
    CCutil_init_timer(&(stat->tot_solve_lp), "tot_solve_lp");
    CCutil_init_timer(&(stat->tot_pricing), "tot_pricing");
    CCutil_init_timer(&(stat->tot_heuristic), "tot_heuristic");
    CCutil_init_timer(&(stat->tot_reduce_cost_fixing),
                      "tot_reduce_cost_fixing");

    stat->mip_nb_vars = 0;
    stat->mip_nb_constr = 0;
    stat->mip_obj_bound = 0.0;
    stat->mip_obj_bound_lp = 0.0;
    stat->mip_rel_gap = 0.0;
    stat->mip_run_time = 110.0;
    stat->mip_status = 0;
    stat->mip_nb_iter_simplex = 0;
    stat->mip_nb_nodes = 0;
    stat->mip_reduced_cost_fixing = 0;

    stat->real_time_build_dd = 0.0;
    stat->real_time_total = getRealTime();
    stat->real_time_branch_and_bound = 0.0;
    stat->real_time_strong_branching = 0.0;
    stat->real_time_lb_root = 0.0;
    stat->real_time_lb = 0.0;
    stat->real_time_pricing = 0.0;
    stat->real_time_heuristic = 0.0;
    CCutil_start_timer(&(stat->tot_cputime));
    stat->first_size_graph = 0;
    stat->size_graph_after_reduced_cost_fixing = 0;
}

void statistics_free(Statistics* stat) {}
