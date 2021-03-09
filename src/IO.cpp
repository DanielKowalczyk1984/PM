#include <boost/date_time/gregorian/gregorian.hpp>
#include "MIP_defs.hpp"
#include "Statistics.h"
#include "wctprivate.h"

void Problem::to_csv() {
    int   val = 0;
    FILE* file = (FILE*)NULL;
    char* file_name = CC_SAFE_MALLOC(128, char);
    auto  d =
        boost::gregorian::date(boost::gregorian::day_clock::universal_day());
    // date;
    stat.real_time_total = getRealTime() - stat.real_time_total;
    CCutil_stop_timer(&(stat.tot_cputime), 0);

    sprintf(file_name, "CG_overall_%d%02d%02d.csv", d.year(), d.month(),
            d.day());

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
            stat.tot_reduce_cost_fixing.name, "rel_error", "global_lower_bound",
            "global_upper_bound", "first_rel_error", "nb_generated_col", "date",
            "nb_iterations_rvnd", "stabilization", "alpha", "pricing_solver",
            "n", "m", "first_size_graph", "size_after_reduced_cost",
            "mip_nb_vars", "mip_nb_constr", "mip_obj_bound", "mip_obj_bound_lp",
            "mip_rel_gap", "mip_run_time", "mip_status", "mip_nb_iter_simplex",
            "mip_nb_nodes");
    }

    // for (int i = MIP_Attr_Run_Time; i <= MIP_Attr_Nb_Nodes &&
    // parms.mip_solver;
    //      i++) {
    //     pd->get_mip_statistics(static_cast<MIP_Attr>(i));
    // }

    fprintf(file,
            "%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%d,%u/"
            "%u/%u,%d,%d,%f,%d,%d,%d,%lu,%lu,%d,%d,%f,%f,%f,%f,%d,%f,%f\n",
            stat.pname.c_str(), stat.real_time_total, stat.tot_cputime.cum_zeit,
            stat.tot_lb.cum_zeit, stat.tot_lb_root.cum_zeit,
            stat.tot_heuristic.cum_zeit, stat.tot_build_dd.cum_zeit,
            stat.tot_pricing.cum_zeit, stat.tot_reduce_cost_fixing.cum_zeit,
            stat.rel_error, stat.global_lower_bound, stat.global_upper_bound,
            stat.root_rel_error, stat.nb_generated_col, d.day(), d.month(),
            d.year(), parms.nb_iterations_rvnd, parms.stab_technique,
            parms.alpha, parms.pricing_solver, nb_jobs, parms.nb_machines,
            stat.first_size_graph, stat.size_graph_after_reduced_cost_fixing,
            stat.mip_nb_vars, stat.mip_nb_constr, stat.mip_obj_bound,
            stat.mip_obj_bound_lp, stat.mip_rel_gap, stat.mip_run_time,
            stat.mip_status, stat.mip_nb_iter_simplex, stat.mip_nb_nodes);
    fclose(file);
CLEAN:
    CC_FREE(file_name, char);
}

int Problem::print_to_screen() {
    int val = 0;

    switch (status) {
        case no_sol:
            fmt::print(
                "We didn't decide if this instance is feasible or "
                "infeasible\n");
            break;

        case feasible:
        case lp_feasible:
        case meta_heuristic:
            fmt::print(
                "A suboptimal schedule with relative error {} is found.\n",
                static_cast<double>(global_upper_bound - global_lower_bound) /
                    (global_lower_bound));
            break;

        case optimal:
            fmt::print("The optimal schedule is found.\n");
            break;
    }

    fmt::print(
        "Compute_schedule took {} seconds(tot_scatter_search {}, "
        "tot_branch_and_bound {}, tot_lb_lp_root {}, tot_lb_lp {}, tot_lb {}, "
        "tot_pricing {}, tot_build_dd {}) and {} seconds in real time\n",
        stat.tot_cputime.cum_zeit, stat.tot_heuristic.cum_zeit,
        stat.tot_branch_and_bound.cum_zeit, stat.tot_lb_root.cum_zeit,
        stat.tot_lb.cum_zeit, stat.tot_lb.cum_zeit, stat.tot_pricing.cum_zeit,
        stat.tot_build_dd.cum_zeit, stat.real_time_total);
    return val;
}
