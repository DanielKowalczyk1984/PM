#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <cstddef>  // for size_t
#include <string>   // for string
#include "util.h"   // for CCutil_timer
struct Parms;

struct Statistics {
    int    global_upper_bound;
    int    global_lower_bound;
    double rel_error;

    int    root_upper_bound;
    int    root_lower_bound;
    double root_rel_error;

    size_t nb_generated_col;
    size_t nb_generated_col_root;
    size_t first_size_graph;
    size_t size_graph_after_reduced_cost_fixing;

    enum TimerType {
        build_dd_timer,
        cputime_timer,
        bb_timer,
        lb_root_timer,
        lb_timer,
        solve_lp_timer,
        pricing_timer,
        heuristic_timer,
        reduced_cost_fixing_timer
    };

    CCutil_timer tot_build_dd;
    CCutil_timer tot_cputime;
    CCutil_timer tot_branch_and_bound;
    CCutil_timer tot_strong_branching;
    CCutil_timer tot_lb_root;
    CCutil_timer tot_lb;
    CCutil_timer tot_solve_lp;
    CCutil_timer tot_pricing;
    CCutil_timer tot_heuristic;
    CCutil_timer tot_reduce_cost_fixing;

    double real_time_build_dd;
    double real_time_total;
    double real_time_branch_and_bound;
    double real_time_strong_branching;
    double real_time_lb_root;
    double real_time_lb;
    double real_time_solve_lp;
    double real_time_pricing;
    double real_time_heuristic;

    int    mip_nb_vars;
    int    mip_nb_constr;
    double mip_obj_bound;
    double mip_obj_bound_lp;
    double mip_rel_gap;
    double mip_run_time;
    int    mip_status;
    double mip_nb_iter_simplex;
    double mip_nb_nodes;
    int    mip_reduced_cost_fixing;

    std::string pname;

    void   start_resume_timer(TimerType _type);
    void   suspend_timer(TimerType _type);
    double total_timer(TimerType _type);

    Statistics(const Parms& parms);
};

#endif  // __STATISTICS_H__