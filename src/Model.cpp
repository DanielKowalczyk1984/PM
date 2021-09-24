#include <fmt/core.h>
#include <OsiGrbSolverInterface.hpp>
#include <array>
#include <cstdio>                                      // for fflush, size_t
#include <memory>                                      // for shared_ptr
#include <range/v3/algorithm/for_each.hpp>             // for for_each, for_...
#include <range/v3/functional/identity.hpp>            // for identity
#include <range/v3/iterator/basic_iterator.hpp>        // for basic_iterator
#include <range/v3/iterator/unreachable_sentinel.hpp>  // for operator==
#include <range/v3/numeric/iota.hpp>                   // for iota, iota_fn
#include <range/v3/range/conversion.hpp>               // for to_vector
#include <range/v3/range_fwd.hpp>                      // for move
#include <range/v3/view/enumerate.hpp>                 // for enumerate_fn
#include <range/v3/view/take.hpp>                      // for take_view, take
#include <range/v3/view/transform.hpp>                 // for transform
#include <range/v3/view/view.hpp>                      // for operator|, vie...
#include <range/v3/view/zip.hpp>                       // for zip_view
#include <range/v3/view/zip_with.hpp>                  // for iter_zip_with_...
#include <span>                                        // for span
#include <vector>                                      // for vector
#include "Column.h"                                    // for ScheduleSet
#include "NodeData.h"                                  // for NodeData
#include "PricerSolverBase.hpp"                        // for PricerSolverBase
#include "Solution.hpp"                                // for Sol
#include "gurobi_c++.h"
#include "gurobi_c.h"  // for GRB_INFINITY
#include "lp.h"        // for lp_interface_g...

int NodeData::add_lhs_column_to_rmp(double cost) {
    id_row.clear();
    coeff_row.clear();

    for (auto const&& [j, it] : lhs_coeff | ranges::views::enumerate) {
        if (std::abs(it) > EPS_BOUND) {
            id_row.emplace_back(j);
            coeff_row.emplace_back(it);
        }
    }

    // lp_interface_addcol(RMP.get(), static_cast<int>(id_row.size()),
    //                     id_row.data(), coeff_row.data(), cost, 0.0,
    //                     GRB_INFINITY, lp_interface_CONT, nullptr);

    osi_rmp->addCol(static_cast<int>(id_row.size()), id_row.data(),
                    coeff_row.data(), 0.0, osi_rmp->getInfinity(), cost);

    return 0;
}

int NodeData::add_lhs_column_to_rmp(double                     cost,
                                    const std::vector<double>& _lhs) {
    id_row.clear();
    coeff_row.clear();

    for (auto const&& [j, it] : _lhs | ranges::views::enumerate) {
        if (std::abs(it) > EPS_BOUND) {
            id_row.emplace_back(j);
            coeff_row.emplace_back(it);
        }
    }

    // lp_interface_addcol(RMP.get(), static_cast<int>(id_row.size()),
    //                     id_row.data(), coeff_row.data(), cost, 0.0,
    //                     GRB_INFINITY, lp_interface_CONT, nullptr);

    osi_rmp->addCol(static_cast<int>(id_row.size()), id_row.data(),
                    coeff_row.data(), 0.0, osi_rmp->getInfinity(), cost);

    return 0;
}

int NodeData::add_column_to_rmp(const Column* set) {
    auto var_ind = 0;
    auto cval = 0.0;
    auto cost = set->total_weighted_completion_time;
    // auto* lp = RMP.get();

    // lp_interface_get_nb_cols(lp, &(nb_cols));
    // var_ind = nb_cols;
    // lp_interface_addcol(lp, 0, nullptr, nullptr, cost, 0.0, GRB_INFINITY,
    //                     lp_interface_CONT, nullptr);

    osi_rmp->addCol(0, nullptr, nullptr, 0.0, osi_rmp->getInfinity(), cost);
    solver->compute_lhs(*set, lhs_coeff.data());
    add_lhs_column_to_rmp(set->total_weighted_completion_time);

    // ranges::for_each(
    //     set->job_list,
    //     [&](auto ind) {
    //         lp_interface_getcoeff(lp, &ind, &var_ind, &cval);
    //         cval += 1.0;
    //         lp_interface_chgcoeff(lp, 1, &ind, &var_ind, &cval);
    //     },
    //     [](const auto tmp) { return static_cast<int>(tmp->job); });

    // auto row_ind = static_cast<int>(nb_jobs);
    // cval = -1.0;
    // lp_interface_chgcoeff(lp, 1, &row_ind, &var_ind, &cval);

    return 0;
}

void NodeData::create_assignment_constraints() {
    id_assignment_constraint = osi_rmp->getNumRows();
    std::vector<int>    start(nb_jobs + 1, 0);
    std::vector<double> rhs_tmp(nb_jobs, 1.0);
    std::vector<double> ub_rhs(nb_jobs, osi_rmp->getInfinity());

    osi_rmp->addRows(static_cast<int>(nb_jobs), start.data(), nullptr, nullptr,
                     rhs_tmp.data(), ub_rhs.data());
}

void NodeData::create_convex_contraint() {
    id_convex_constraint = osi_rmp->getNumRows();
    auto               rhs_lb = -static_cast<double>(nb_machines);
    auto               rhs_ub = osi_rmp->getInfinity();
    std::array<int, 2> start = {0, 0};

    osi_rmp->addRows(1, start.data(), nullptr, nullptr, &rhs_lb, &rhs_ub);

    id_valid_cuts = osi_rmp->getNumRows();
    nb_rows = id_valid_cuts;
}

void NodeData::create_artificial_cols() {
    id_art_var_assignment = osi_rmp->getNumCols();
    id_art_var_convex = nb_jobs;
    id_art_var_cuts = nb_jobs + 1;
    id_next_var_cuts = id_art_var_cuts;
    auto nb_vars = nb_jobs + 1 + max_nb_cuts;
    id_pseudo_schedules = static_cast<int>(nb_vars);

    std::vector<double> lb(nb_vars, 0.0);
    std::vector<double> ub(nb_vars, GRB_INFINITY);
    std::vector<double> obj(nb_vars, 100.0 * (opt_sol.tw + 1));
    std::vector<int> start_vars(nb_vars + 1, static_cast<int>(nb_jobs + 1UL));

    auto nz = static_cast<int>(nb_jobs + 1UL);

    std::vector<int> rows_ind(nz);
    ranges::iota(rows_ind, 0);
    std::vector<double> coeff_vals(nz, 1.0);
    coeff_vals[nb_jobs] = -1.0;
    ranges::iota(start_vars | ranges::views::take(nz), 0);

    osi_rmp->addCols(static_cast<int>(nb_vars), start_vars.data(),
                     rows_ind.data(), coeff_vals.data(), lb.data(), ub.data(),
                     obj.data());

    id_pseudo_schedules = osi_rmp->getNumCols();
}

void NodeData::add_cols_local_pool() {
    // lp_interface_get_nb_rows(RMP.get(), &nb_rows);
    nb_rows = osi_rmp->getNumRows();
    std::vector<double> _lhs(nb_rows);
    ranges::for_each(localColPool, [&](auto& it) {
        solver->compute_lhs(*it.get(), _lhs.data());
        add_lhs_column_to_rmp(it->total_weighted_completion_time, _lhs);
    });
}

int NodeData::build_rmp() {
    /**
     * Set up messegahandler osi solver
     */
    osi_rmp->messageHandler()->setLogLevel(0);
    osi_rmp->messageHandler()->setPrefix(false);
    osi_rmp->setHintParam(OsiDoReducePrint, true, OsiHintTry);

    // std::vector<int>    start(nb_jobs + 1, 0);
    // std::vector<double> rhs_tmp(nb_jobs, 1.0);
    // std::vector<char>   sense(nb_jobs, GRB_GREATER_EQUAL);

    /**
     * add assignment constraints
     */
    // lp_interface_get_nb_rows(RMP.get(), &(id_assignment_constraint));
    // lp_interface_addrows(RMP.get(), static_cast<int>(nb_jobs), 0,
    // start.data(),
    //                      nullptr, nullptr, sense.data(), rhs_tmp.data(),
    //                      nullptr);

    /**
     * add number of machines constraint (convexification)
     */
    // lp_interface_get_nb_rows(RMP.get(), &(id_convex_constraint));
    // lp_interface_addrow(RMP.get(), 0, nullptr, nullptr,
    //                     lp_interface_GREATER_EQUAL,
    //                     -static_cast<double>(nb_machines), nullptr);
    // lp_interface_get_nb_rows(RMP.get(), &(id_valid_cuts));
    // nb_rows = id_valid_cuts;

    /**
     * construct artificial variables in RMP
     */
    // lp_interface_get_nb_cols(RMP.get(), &(id_art_var_assignment));
    // id_art_var_convex = nb_jobs;
    // id_art_var_cuts = nb_jobs + 1;
    // id_next_var_cuts = id_art_var_cuts;
    // auto nb_vars = nb_jobs + 1 + max_nb_cuts;
    // id_pseudo_schedules = static_cast<int>(nb_vars);

    // std::vector<double> lb(nb_vars, 0.0);
    // std::vector<double> ub(nb_vars, GRB_INFINITY);
    // std::vector<double> obj(nb_vars, 100.0 * (opt_sol.tw + 1));
    // std::vector<char>   vtype(nb_vars, GRB_CONTINUOUS);
    // std::vector<int> start_vars(nb_vars + 1, static_cast<int>(nb_jobs +
    // 1UL));

    // auto nz = static_cast<int>(nb_jobs + 1UL);

    // std::vector<int> rows_ind(nz);
    // ranges::iota(rows_ind, 0);
    // std::vector<double> coeff_vals(nz, 1.0);
    // coeff_vals[nb_jobs] = -1.0;
    // ranges::iota(start_vars | ranges::views::take(nz), 0);

    // lp_interface_addcols(RMP.get(), static_cast<int>(nb_vars), nz,
    //                      start_vars.data(), rows_ind.data(),
    //                      coeff_vals.data(), obj.data(), lb.data(), ub.data(),
    //                      vtype.data(), nullptr);
    // lp_interface_get_nb_cols(RMP.get(), &(id_pseudo_schedules));

    create_assignment_constraints();
    create_convex_contraint();
    create_artificial_cols();

    prune_duplicated_sets();
    add_cols_local_pool();

    /**
     * Some aux variables for column generation
     */
    // osi_rmp->initialSolve();
    // pi.resize(nb_jobs + 1, 0.0);
    slack.resize(nb_jobs + 1, 0.0);
    // rhs.resize(nb_jobs + 1, 0.0);
    rhs = std::span<const double>(osi_rmp->getRightHandSide(),
                                  static_cast<size_t>(osi_rmp->getNumRows()));
    // lp_interface_get_rhs(RMP.get(), rhs.data());
    lhs_coeff.resize(nb_jobs + 1, 0.0);
    id_row.reserve(nb_jobs + 1);
    coeff_row.reserve(nb_jobs + 1);

    return 0;
}
