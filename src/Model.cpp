#include <fmt/core.h>
#include <cassert>                                     // for assert
#include <cstdio>                                      // for fflush, size_t
#include <memory>                                      // for shared_ptr
#include <range/v3/algorithm/for_each.hpp>             // for for_each, for_...
#include <range/v3/functional/identity.hpp>            // for identity
#include <range/v3/iterator/basic_iterator.hpp>        // for basic_iterator
#include <range/v3/iterator/unreachable_sentinel.hpp>  // for operator==
#include <range/v3/numeric/iota.hpp>                   // for iota, iota_fn
#include <range/v3/range_fwd.hpp>                      // for move
#include <range/v3/view/drop.hpp>                      // for drop_view, drop
#include <range/v3/view/enumerate.hpp>                 // for enumerate_fn
#include <range/v3/view/take.hpp>                      // for take_view, take
#include <range/v3/view/view.hpp>                      // for operator|, vie...
#include <range/v3/view/zip.hpp>                       // for zip_view
#include <range/v3/view/zip_with.hpp>                  // for iter_zip_with_...
#include <vector>                                      // for vector
#include "Solution.hpp"                                // for Sol
#include "gurobi_c.h"                                  // for GRB_INFINITY
#include "lp.h"                                        // for lp_interface_g...
#include "scheduleset.h"                               // for ScheduleSet
#include "wctprivate.h"                                // for NodeData, EPS_...

int NodeData::add_lhs_scheduleset_to_rmp(Column* set) {
    id_row.clear();
    coeff_row.clear();

    for (auto const&& [j, it] : lhs_coeff | ranges::views::enumerate) {
        if (std::abs(it) > EPS_BOUND) {
            id_row.emplace_back(j);
            coeff_row.emplace_back(it);
        }
    }

    auto cost = set->total_weighted_completion_time;
    lp_interface_addcol(RMP.get(), static_cast<int>(id_row.size()),
                        id_row.data(), coeff_row.data(), cost, 0.0,
                        GRB_INFINITY, lp_interface_CONT, nullptr);

    return 0;
}

int NodeData::add_scheduleset_to_rmp(Column* set) {
    auto  var_ind = 0;
    auto  cval = 0.0;
    auto  cost = set->total_weighted_completion_time;
    auto* lp = RMP.get();

    lp_interface_get_nb_cols(lp, &(nb_cols));
    var_ind = nb_cols;
    lp_interface_addcol(lp, 0, nullptr, nullptr, cost, 0.0, GRB_INFINITY,
                        lp_interface_CONT, nullptr);

    ranges::for_each(
        set->job_list,
        [&](auto ind) {
            lp_interface_getcoeff(lp, &ind, &var_ind, &cval);
            cval += 1.0;
            lp_interface_chgcoeff(lp, 1, &ind, &var_ind, &cval);
        },
        [](const auto tmp) { return static_cast<int>(tmp->job); });

    auto row_ind = static_cast<int>(nb_jobs);
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
    lp_interface_addrows(RMP.get(), static_cast<int>(nb_jobs), 0, start.data(),
                         nullptr, nullptr, sense.data(), rhs_tmp.data(),
                         nullptr);

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
    id_pseudo_schedules = static_cast<int>(nb_vars);

    std::vector<double> lb(nb_vars, 0.0);
    std::vector<double> ub(nb_vars, GRB_INFINITY);
    std::vector<double> obj(nb_vars, 100.0 * (opt_sol.tw + 1));
    std::vector<char>   vtype(nb_vars, GRB_CONTINUOUS);
    std::vector<int> start_vars(nb_vars + 1, static_cast<int>(nb_jobs + 1UL));

    start_vars[nb_vars] = static_cast<int>(nb_jobs + 1);

    std::vector<int> rows_ind(nb_jobs + 1);
    ranges::iota(rows_ind, 0);
    std::vector<double> coeff_vals(nb_jobs + 1, 1.0);
    coeff_vals[nb_jobs] = -1.0;
    ranges::iota(start_vars | ranges::views::take(nb_jobs + 1), 0);

    lp_interface_addcols(RMP.get(), static_cast<int>(nb_vars),
                         static_cast<int>(nb_jobs + 1), start_vars.data(),
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
