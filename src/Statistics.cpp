#include "Statistics.h"

Statistics::Statistics()
    : global_upper_bound(INT_MAX),
      global_lower_bound(0),
      rel_error(DBL_MAX),
      root_upper_bound(INT_MAX),
      root_lower_bound(0.0),
      root_rel_error(DBL_MAX),
      nb_explored_nodes(0),
      nb_generated_nodes(0),
      nb_generated_col(0),
      nb_generated_col_root(0),
      first_size_graph(0),
      size_graph_after_reduced_cost_fixing(0),
      real_time_build_dd(0.0),
      real_time_total(getRealTime()),
      real_time_branch_and_bound(0.0),
      real_time_strong_branching(0.0),
      real_time_lb_root(0.0),
      real_time_lb(0.0),
      real_time_solve_lp(0.0),
      real_time_pricing(0.0),
      real_time_heuristic(0.0),
      mip_nb_vars(0),
      mip_nb_constr(0),
      mip_obj_bound(0.0),
      mip_obj_bound_lp(0.0),
      mip_rel_gap(0.0),
      mip_run_time(110.0),
      mip_status(0),
      mip_nb_iter_simplex(),
      mip_nb_nodes(0),
      mip_reduced_cost_fixing(),
      pname() {
    CCutil_init_timer(&(tot_build_dd), "tot_build_dd");
    CCutil_init_timer(&(tot_cputime), "tot_cputime");
    CCutil_init_timer(&(tot_branch_and_bound), "tot_branch_and_bound");
    CCutil_init_timer(&(tot_strong_branching), "tot_strong_branching");
    CCutil_init_timer(&(tot_lb_root), "tot_lb_root");
    CCutil_init_timer(&(tot_lb), "tot_lb");
    CCutil_init_timer(&(tot_solve_lp), "tot_solve_lp");
    CCutil_init_timer(&(tot_pricing), "tot_pricing");
    CCutil_init_timer(&(tot_heuristic), "tot_heuristic");
    CCutil_init_timer(&(tot_reduce_cost_fixing), "tot_reduce_cost_fixing");
    CCutil_start_timer(&(tot_cputime));
}
