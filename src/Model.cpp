#include <fmt/core.h>
#include <numeric>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/numeric/iota.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/take.hpp>
#include <vector>
#include "Instance.h"
#include "Statistics.h"
#include "lp.h"
#include "scheduleset.h"
#include "wctprivate.h"

int NodeData::grab_integer_solution(std::vector<double> const& x,
                                    double                     tolerance) {
    auto incumbent = 0.0;
    auto tot_weighted = 0;

    lp_interface_objval(RMP.get(), &incumbent);
    lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    assert(nb_cols - id_pseudo_schedules == localColPool.size());

    best_schedule.clear();

    for (const auto&& [i, tmp_schedule] :
         localColPool | ranges::views::enumerate |
             ranges::views::drop(id_pseudo_schedules)) {
        if (x[i] >= 1.0 - tolerance) {
            auto aux = std::make_shared<ScheduleSet>(*tmp_schedule);
            best_schedule.emplace_back(aux);

            tot_weighted += aux->total_weighted_completion_time;

            if (best_schedule.size() > nb_machines) {
                fmt::print(
                    "ERROR: \"Integral\" solution turned out to be not "
                    "integral!\n");
                fflush(stdout);
            }
        }
    }

    /** Write a check function */
    fmt::print("Intermediate schedule:\n");
    fmt::print("with total weight {}\n", tot_weighted);
    assert(fabs((double)tot_weighted - incumbent) <= EPS);

    if (tot_weighted < upper_bound) {
        upper_bound = tot_weighted;
        best_objective = tot_weighted;
    }

    if (upper_bound == lower_bound) {
        status = finished;
    }

    return 0;
}

int NodeData::add_lhs_scheduleset_to_rmp(ScheduleSet* set) {
    id_row.clear();
    coeff_row.clear();

    for (auto const&& [j, it] : lhs_coeff | ranges::views::enumerate) {
        if (std::abs(it) > EPS_BOUND) {
            id_row.emplace_back(j);
            coeff_row.emplace_back(it);
        }
    }

    int len = id_row.size();

    auto cost = static_cast<double>(set->total_weighted_completion_time);
    lp_interface_addcol(RMP.get(), len, id_row.data(), coeff_row.data(), cost,
                        0.0, GRB_INFINITY, lp_interface_CONT, nullptr);

    return 0;
}

int NodeData::add_scheduleset_to_rmp(ScheduleSet* set) {
    auto  row_ind = 0;
    auto  var_ind = 0;
    auto  cval = 0.0;
    auto  cost = static_cast<double>(set->total_weighted_completion_time);
    auto* lp = RMP.get();

    lp_interface_get_nb_cols(lp, &(nb_cols));
    var_ind = nb_cols;
    lp_interface_addcol(lp, 0, nullptr, nullptr, cost, 0.0, GRB_INFINITY,
                        lp_interface_CONT, nullptr);

    for (auto& it : set->job_list) {
        row_ind = it->job;
        lp_interface_getcoeff(lp, &row_ind, &var_ind, &cval);
        cval += 1.0;
        lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);
    }

    row_ind = nb_jobs;
    cval = -1.0;
    lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);

    return 0;
}

int NodeData::build_rmp() {
    std::vector<int>    start(nb_jobs, 0);
    std::vector<double> rhs_tmp(nb_jobs, 1.0);
    std::vector<char>   sense(nb_jobs, GRB_GREATER_EQUAL);

    /**
     * add assignment constraints
     */
    lp_interface_get_nb_rows(RMP.get(), &(id_assignment_constraint));
    lp_interface_addrows(RMP.get(), nb_jobs, 0, start.data(), nullptr, nullptr,
                         sense.data(), rhs_tmp.data(), nullptr);

    /**
     * add number of machines constraint (convexification)
     */
    lp_interface_get_nb_rows(RMP.get(), &(id_convex_constraint));
    lp_interface_addrow(RMP.get(), 0, nullptr, nullptr,
                        lp_interface_GREATER_EQUAL,
                        -static_cast<double>(nb_machines), nullptr);
    lp_interface_get_nb_rows(RMP.get(), &(id_valid_cuts));
    nb_rows = id_valid_cuts;

    /**
     * construct artificial variables in RMP
     */
    lp_interface_get_nb_cols(RMP.get(), &(id_art_var_assignment));
    id_art_var_convex = nb_jobs;
    id_art_var_cuts = nb_jobs + 1;
    id_next_var_cuts = id_art_var_cuts;
    auto nb_vars = nb_jobs + 1 + max_nb_cuts;
    id_pseudo_schedules = nb_vars;

    std::vector<double> lb(nb_vars, 0.0);
    std::vector<double> ub(nb_vars, GRB_INFINITY);
    std::vector<double> obj(nb_vars, 100.0 * (opt_sol.tw + 1));
    std::vector<char>   vtype(nb_vars, GRB_CONTINUOUS);
    std::vector<int>    start_vars(nb_vars + 1, nb_jobs + 1UL);

    start_vars[nb_vars] = nb_jobs + 1;

    std::vector<int> rows_ind(nb_jobs + 1);
    ranges::iota(rows_ind, 0);
    std::vector<double> coeff_vals(nb_jobs + 1, 1.0);
    coeff_vals[nb_jobs] = -1.0;
    ranges::iota(start_vars.begin(), start_vars.begin() + nb_jobs + 1, 0);

    lp_interface_addcols(RMP.get(), nb_vars, nb_jobs + 1, start_vars.data(),
                         rows_ind.data(), coeff_vals.data(), obj.data(),
                         lb.data(), ub.data(), vtype.data(), nullptr);
    lp_interface_get_nb_cols(RMP.get(), &(id_pseudo_schedules));

    /** add columns from localColPool */
    prune_duplicated_sets();
    ranges::for_each(localColPool,
                     [&](auto& it) { add_scheduleset_to_rmp(it.get()); });

    /**
     * Some aux variables for column generation
     */

    pi.resize(nb_jobs + 1, 0.0);
    slack.resize(nb_jobs + 1, 0.0);
    rhs.resize(nb_jobs + 1, 0.0);
    lp_interface_get_rhs(RMP.get(), rhs.data());
    lhs_coeff.resize(nb_jobs + 1, 0.0);
    id_row.reserve(nb_jobs + 1);
    coeff_row.reserve(nb_jobs + 1);

    return 0;
}
