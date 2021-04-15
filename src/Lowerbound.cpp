#include <bits/ranges_algo.h>
#include <fmt/core.h>
#include <range/v3/action/remove_if.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>
#include <vector>
#include "Instance.h"
#include "Job.h"
#include "PricerSolverBase.hpp"
#include "Statistics.h"
#include "gurobi_c.h"
#include "lp.h"
#include "scheduleset.h"
#include "util.h"
#include "wctprivate.h"

/** Help function for column generation */
void NodeData::print_ages() {
    fmt::print("AGES:");

    std::ranges::for_each(localColPool,
                          [](auto const& it) { fmt::print(" {}", it->age); });

    fmt::print("\n");
}

int NodeData::grow_ages() {
    int val = 0;
    lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    assert(((nb_cols - id_pseudo_schedules) == localColPool.size()));
    if (!localColPool.empty()) {
        column_status.resize(localColPool.size());
        val = lp_interface_basis_cols(RMP.get(), column_status.data(),
                                      id_pseudo_schedules);
        // CCcheck_val_2(val, "Failed in lp_interface_basis_cols");
        zero_count = 0;

        std::ranges::for_each(
            localColPool | ranges::views::enumerate, [&](auto&& it) {
                if (column_status[it.first] == lp_interface_LOWER ||
                    column_status[it.first] == lp_interface_FREE) {
                    it.second->age++;

                    if (it.second->age > retirementage) {
                        zero_count++;
                    }
                } else {
                    it.second->age = 0;
                }
            });
    }

    return val;
}

int NodeData::delete_unused_rows() {
    int              val = 0;
    std::vector<int> del_indices{};

    lp_interface_get_nb_rows(RMP.get(), &nb_rows);
    // assert(nb_rows == pd->slack->len);
    lp_interface_slack(RMP.get(), slack.data());

    int it = id_valid_cuts;
    int first_del = -1;
    int last_del = -1;
    for (int i = id_valid_cuts; i < nb_rows; i++) {
        if (std::fabs(slack[i]) < EPS) {
            if (first_del != -1) {
                val = delete_unused_rows_range(first_del, last_del);
                it = it - (last_del - first_del);
                first_del = last_del = -1;
            } else {
                it++;
            }
        } else {
            if (first_del == -1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            it++;
        }
    }

    if (first_del != -1) {
        delete_unused_rows_range(first_del, last_del);
    }

    call_update_rows_coeff();

    return val;
}

int NodeData::delete_old_schedules() {
    int val = 0;
    int min_numdel = floor(nb_jobs * min_nb_del_row_ratio);
    /** pd->zero_count can be deprecated! */
    zero_count = 0;

    std::ranges::for_each(localColPool, [&](auto const& it) {
        if (it->age > 0) {
            zero_count++;
        }
    });

    if (zero_count > min_numdel) {
        int              iter = 0;
        std::vector<int> dellist{};
        lp_interface_get_nb_cols(RMP.get(), &nb_cols);
        assert(nb_cols - id_pseudo_schedules == localColPool.size());

        localColPool |= ranges::actions::remove_if([&](const auto& it) {
            auto val = false;
            if (it->age > retirementage) {
                dellist.emplace_back(iter + id_pseudo_schedules);
                val = true;
            }
            iter++;
            return val;
        });

        lp_interface_delete_cols_array(RMP.get(), dellist.data(),
                                       dellist.size());

        if (dbg_lvl() > 1) {
            fmt::print("Deleted {} out of {} columns with age > {}.\n",
                       zero_count, localColPool.size(), retirementage);
        }
        lp_interface_get_nb_cols(RMP.get(), &nb_cols);
        assert(localColPool.size() == nb_cols - id_pseudo_schedules);
        zero_count = 0;
    }

    return val;
}

int NodeData::delete_infeasible_schedules() {
    auto count = localColPool.size();
    /** pd->zero_count can be deprecated! */
    zero_count = 0;

    int iter = 0;
    lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    assert(nb_cols - id_pseudo_schedules == count);
    std::vector<int> dellist{};

    // std::erase_if(localColPool, );

    localColPool |= ranges::actions::remove_if([&](auto& it) {
        auto val = false;
        if (!it->del) {
            dellist.emplace_back(iter + id_pseudo_schedules);
            val = true;
        }
        iter++;
        return val;
    });

    lp_interface_delete_cols_array(RMP.get(), dellist.data(), dellist.size());

    if (dbg_lvl() > 1) {
        fmt::print(
            "Deleted {} out of {} columns(infeasible columns after "
            "reduce "
            "cost "
            "fixing).\n",
            dellist.size(), count);
    }

    lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    assert(localColPool.size() == nb_cols - id_pseudo_schedules);
    if (dbg_lvl() > 1) {
        fmt::print("number of cols = {}\n", nb_cols - id_pseudo_schedules);
    }

    if (!dellist.empty()) {
        solve_relaxation();
    }

    return 0;
}

// void g_make_pi_feasible(gpointer data, gpointer user_data) {
//     ScheduleSet* x = static_cast<ScheduleSet*>(data);
//     NodeData*    pd = static_cast<NodeData*>(user_data);
//     Job*         tmp_j = nullptr;

//     double colsum = .0;

//     for (guint i = 0; i < x->job_list->len; ++i) {
//         tmp_j = static_cast<Job*>(g_ptr_array_index(x->job_list, i));
//         if (signbit(pd->pi[tmp_j->job])) {
//             pd->pi[tmp_j->job] = 0.0;
//         }

//         colsum += pd->pi[tmp_j->job];
//         colsum = nextafter(colsum, DBL_MAX);
//     }

//     if (!signbit(pd->pi[pd->nb_jobs])) {
//         pd->pi[pd->nb_jobs] = 0.0;
//     }

//     colsum += pd->pi[pd->nb_jobs];
//     colsum = nextafter(colsum, DBL_MAX);

//     if (colsum > x->total_weighted_completion_time) {
//         double newcolsum = .0;
//         for (guint i = 0; i < x->job_list->len; ++i) {
//             tmp_j = static_cast<Job*>(g_ptr_array_index(x->job_list, i));
//             pd->pi[tmp_j->job] /= colsum;
//             pd->pi[tmp_j->job] *= x->total_weighted_completion_time;
//             newcolsum += pd->pi[tmp_j->job];
//         }

//         pd->pi[pd->nb_jobs] /= colsum;
//         pd->pi[pd->nb_jobs] *= x->total_weighted_completion_time;
//         newcolsum += pd->pi[pd->nb_jobs];

//         if (dbg_lvl() > 1) {
//             fmt::print(
//                 R"(Decreased column sum of {} from  {:30.20f} to
//                 {:30.20f}
// )",
//                 x->id, colsum, newcolsum);
//         }
//     }
// }

// MAYBE_UNUSED
// void make_pi_feasible(NodeData* pd) {
//     g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible, pd);

//     std::for_each()
// }

// void g_make_pi_feasible_farkas(gpointer data, gpointer user_data) {
//     ScheduleSet* x = static_cast<ScheduleSet*>(data);
//     NodeData*    pd = static_cast<NodeData*>(user_data);

//     double colsum = .0;

//     for (guint i = 0; i < x->job_list->len; ++i) {
//         double tmp = pd->pi[i];
//         if (signbit(tmp)) {
//             tmp = 0.0;
//         }

//         colsum += tmp;
//         colsum = nextafter(colsum, DBL_MAX);
//     }

//     colsum += pd->pi[pd->nb_jobs];

//     if (colsum > x->total_weighted_completion_time) {
//         double newcolsum = .0;
//         for (guint i = 0; i < x->job_list->len; ++i) {
//             double tmp = pd->pi[i];
//             tmp /= colsum;
//             tmp *= x->total_weighted_completion_time;
//             newcolsum += tmp;
//         }

//         double tmp = pd->pi[pd->nb_jobs];
//         tmp /= colsum;
//         tmp *= x->total_weighted_completion_time;
//         newcolsum += tmp;

//         if (dbg_lvl() > 1) {
//             fmt::print(
//                 R"(Decreased column sum of {} from  {:30.20f} to
//                 {:30.20f}
// )",
//                 x->id, colsum, newcolsum);
//         }
//     }
// }

// MAYBE_UNUSED
// void NodeData::make_pi_feasible_farkas_pricing() {
//     g_ptr_array_foreach(localColPool, g_make_pi_feasible_farkas, this);
// }

int NodeData::compute_objective() {
    LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    assert(nb_rows == pi.size());

    LP_lower_bound_dual = ranges::inner_product(pi, rhs, 0.0);
    LP_lower_bound_dual -= EPS_BOUND;

    /** Get the LP lower bound and compute the lower bound of WCT */
    lp_interface_objval(RMP.get(), &(LP_lower_bound));
    LP_lower_bound -= EPS_BOUND;
    // CCcheck_val_2(val, "lp_interface_objval failed");
    lower_bound = (ceil(LP_lower_bound_dual) < ceil(LP_lower_bound))
                      ? static_cast<int>(ceil(LP_lower_bound_dual))
                      : static_cast<int>(ceil(LP_lower_bound));
    LP_lower_bound_BB = std::min(LP_lower_bound, LP_lower_bound_dual);
    LP_lower_min = std::min(LP_lower_min, LP_lower_bound_BB);

    if (iterations % (nb_jobs) == 0 && dbg_lvl() > 0) {
        fmt::print(
            "Current primal LP objective: {:19.16f}  (LP_dual-bound "
            "{:19.16f}, "
            "lowerbound = {}, eta_in = {}, eta_out = {}).\n",
            LP_lower_bound + instance.off, LP_lower_bound_dual + instance.off,
            lower_bound + instance.off,
            solver_stab->get_eta_in() + instance.off,
            LP_lower_bound + instance.off);
    }

    return 0;
}

int NodeData::solve_relaxation() {
    int    val = 0;
    int    status = 0;
    double real_time_solve_lp = 0.0;

    /** Compute LP relaxation */
    real_time_solve_lp = getRealTime();
    stat.start_resume_timer(Statistics::solve_lp_timer);
    val = lp_interface_optimize(RMP.get(), &status);
    // CCcheck_val_2(val, "lp_interface_optimize failed");
    stat.suspend_timer(Statistics::solve_lp_timer);
    real_time_solve_lp = getRealTime() - real_time_solve_lp;
    stat.real_time_solve_lp += real_time_solve_lp;

    if (dbg_lvl() > 1) {
        fmt::print("Simplex took {} seconds.\n", real_time_solve_lp);
        fflush(stdout);
    }

    if (dbg_lvl() > 1) {
        print_ages();
    }

    switch (status) {
        case LP_INTERFACE_OPTIMAL:
            /** grow ages of the different columns */
            grow_ages();
            /** get the dual variables and make them feasible */
            lp_interface_pi(RMP.get(), pi.data());
            /** Compute the objective function */
            compute_objective();
            break;

        case LP_INTERFACE_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = lp_interface_pi_inf(RMP.get(), pi.data());
            break;
    }

CLEAN:

    return val;
}

int NodeData::compute_lower_bound() {
    auto j = 0;
    auto val = 0;
    auto has_cols = 1;
    auto has_cuts = 0;
    auto nb_non_improvements = 0;
    auto status_RMP = GRB_LOADED;
    auto real_time_pricing = 0.0;

    if (dbg_lvl() > 1) {
        fmt::print(
            R"(Starting compute_lower_bound with lb {} and ub %d at depth {}
)",
            lower_bound, upper_bound, depth);
    }

    stat.start_resume_timer(Statistics::lb_timer);

    /**
     * Construction of new solution if localPoolColPool is empty
     */
    if (localColPool.empty()) {
        add_solution_to_colpool(opt_sol);
    }

    if (!RMP) {
        val = build_rmp();
    }

    check_schedules();
    delete_infeasible_schedules();

    // solve_relaxation(problem, pd);
    do {
        has_cols = 1;
        has_cuts = 0;
        while ((iterations < NB_CG_ITERATIONS) && has_cols &&
               stat.total_timer(Statistics::cputime_timer) <=
                   parms.branching_cpu_limit) {
            /**
             * Delete old columns
             */
            if (zero_count > nb_jobs * min_nb_del_row_ratio &&
                status == GRB_OPTIMAL) {
                delete_old_schedules();
            }
            solve_relaxation();

            /**
             * Solve the pricing problem
             */
            real_time_pricing = getRealTime();
            stat.start_resume_timer(Statistics::pricing_timer);
            val = lp_interface_status(RMP.get(), &status_RMP);

            switch (status_RMP) {
                case GRB_OPTIMAL:
                    iterations++;
                    status = infeasible;

                    val = solve_pricing();
                    break;

                case GRB_INFEASIBLE:
                    solve_farkas_dbl();
                    break;
            }

            stat.suspend_timer(Statistics::pricing_timer);
            real_time_pricing = getRealTime() - real_time_pricing;
            stat.real_time_pricing += real_time_pricing;

            switch (status_RMP) {
                case GRB_OPTIMAL:
                    has_cols = (solver_stab->stopping_criteria() &&
                                solver_stab->get_eta_in() <
                                    upper_bound - 1.0 + EPS_BOUND);
                    // nb_new_sets = 0;
                    // || nb_non_improvements > 5;  // ||
                    // (ceil(eta_in - 0.00001) >= eta_out);

                    break;

                case GRB_INFEASIBLE:
                    has_cols = (nb_new_sets == 0);
                    // nb_new_sets = 0;
                    break;
            }
        }

        switch (status_RMP) {
            case GRB_OPTIMAL:
                if (dbg_lvl() > 1) {
                    fmt::print(
                        R"(Found lb = {} ({}) upper_bound = {} (iterations = {}).
)",
                        lower_bound, LP_lower_bound, upper_bound, iterations);
                }

                /**
                 * Compute the objective function
                 */
                // retirementage = 0;
                // delete_old_schedules(pd);
                solve_relaxation();
                // double obj;
                // lp_interface_objval(pd->RMP, &obj);
                // printf("test objval = %f\n", obj);
                // check_schedules(pd);
                // delete_infeasible_schedules(pd);
                // solve_relaxation(pd);
                // lp_interface_objval(pd->RMP, &obj);
                // printf("test objval = %f\n", obj);
                // printf("----------------\n");
                // compute_objective(pd);
                if (!localColPool.empty()) {
                    construct_lp_sol_from_rmp();
                    // CCcheck_val_2(val, "Failed in construct lp sol
                    // from rmp\n"); solve_relaxation(pd);
                    // delete_old_schedules(pd); delete_unused_rows(pd);
                    // solve_relaxation(pd);
                    // construct_lp_sol_from_rmp(pd);
                    // if (!call_is_integer_solution(pd->solver)) {
                    // has_cuts = (generate_cuts(pd) > 0);
                    // has_cuts = 0;
                    // call_update_duals(pd->solver_stab);
                    // lp_interface_write(pd->RMP, "test.lp");
                    // }
                }
                break;

            case GRB_INFEASIBLE:
                status = infeasible;
                lp_interface_write(RMP.get(), "infeasible_RMP.lp");
                lp_interface_compute_IIS(RMP.get());
        }
    } while (false);

    if (iterations < NB_CG_ITERATIONS &&
        stat.total_timer(Statistics::cputime_timer) <=
            parms.branching_cpu_limit) {
    } else {
        switch (status_RMP) {
            case GRB_OPTIMAL:
                status = LP_bound_computed;
                break;

            case GRB_INFEASIBLE:
                status = infeasible;
                break;
        }
    }
    // } while (depth == 1);

    if (depth == 0) {
        stat.global_lower_bound =
            std::max(lower_bound + instance.off, stat.global_lower_bound);
        stat.root_lower_bound = stat.global_lower_bound;
        stat.root_upper_bound = stat.global_upper_bound;
        stat.root_rel_error = static_cast<double>(stat.global_upper_bound -
                                                  stat.global_lower_bound) /
                              (stat.global_lower_bound + EPS);
        stat.size_graph_after_reduced_cost_fixing = solver->get_nb_vertices();
        stat.nb_generated_col_root = iterations;
    }

    stat.nb_generated_col += iterations;
    stat.suspend_timer(Statistics::lb_timer);

CLEAN:
    return val;
}

int NodeData::print_x() {
    int val = 0;
    int nb_cols = 0;
    int status = 0;

    val = lp_interface_status(RMP.get(), &status);
    // CCcheck_val_2(val, "Failed in lp_interface_status");

    switch (status) {
        case GRB_OPTIMAL:
            val = lp_interface_get_nb_cols(RMP.get(), &nb_cols);
            assert(localColPool.size() == nb_cols - id_pseudo_schedules);
            lambda.resize(nb_cols - id_pseudo_schedules, 0.0);
            val = lp_interface_x(RMP.get(), lambda.data(), id_pseudo_schedules);

            for (auto i = 0UL; auto& it : localColPool) {
                if (lambda[i] > EPS) {
                    std::ranges::for_each(it->job_list, [](const Job* j) {
                        fmt::print("{} ", j->job);
                    });
                    fmt::print("\n");
                }
                ++i;
            }
            break;
    }

CLEAN:

    return val;
}

int NodeData::check_schedules() {
    int status = 0;

    lp_interface_status(RMP.get(), &status);
    lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    assert(nb_cols - id_pseudo_schedules == localColPool.size());
    if (dbg_lvl() > 1) {
        fmt::print("number of cols check {}\n", nb_cols - id_pseudo_schedules);
    }
    for (auto& col : localColPool) {
        if (check_schedule_set(col.get())) {
            col->del = 1;
        } else {
            col->del = 0;
        }
    }

    return 0;
}