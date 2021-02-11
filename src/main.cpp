////////////////////////////////////////////////////////////////
//                                                            //
//  main.c                                                    //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <fmt/core.h>
#include "branch-and-bound/btree.h"
#include "parms.h"
#include "solver.h"
#include "util.h"
#include "wct.h"
#include "wctprivate.h"

int main(int ac, const char** av) {
    int     val = 0;
    Problem problem(ac, av);

    NodeData*   root = problem.root_pd;
    Parms*      parms = &(problem.parms);
    Statistics* statistics = &(problem.stat);

    /**
     *@brief Finding heuristic solutions to the problem or start without
     *feasible solutions
     *
     */
    if (parms->use_heuristic) {
        heuristic(&problem);
    } else {
        problem.opt_sol =
            solution_alloc(problem.intervals->len, root->nb_machines,
                           root->nb_jobs, problem.off);
        Solution* sol = problem.opt_sol;
        // CCcheck_NULL_2(sol, "Failed to allocate memory");
        val = construct_edd(&problem, sol);
        // CCcheck_val_2(val, "Failed construct edd");
        fmt::print("Solution Constructed with EDD heuristic:\n");
        solution_print(sol);
        solution_canonical_order(sol, problem.intervals);
        fmt::print("Solution in canonical order: \n");
        solution_print(sol);
    }

    /**
     * @brief Build DD at the root node
     *
     */
    if (parms->pricing_solver < dp_solver) {
        CCutil_start_timer(&(statistics->tot_build_dd));
        root->solver =
            newSolver(root->jobarray, root->nb_machines, root->ordered_jobs,
                      parms, root->H_max, NULL, problem.opt_sol->tw);
        CCutil_stop_timer(&(statistics->tot_build_dd), 0);
    } else {
        root->solver = newSolverDp(root->jobarray, root->nb_machines,
                                   root->H_max, parms, problem.opt_sol->tw);
    }
    statistics->first_size_graph = get_nb_vertices(root->solver);

    /**
     * @brief Initial stabilization
     *
     */
    root->solver_stab = new_pricing_stabilization(root->solver, parms);
    root->stat = &(problem.stat);
    root->opt_sol = problem.opt_sol;

    /**
     * @brief Solve initial relaxation
     *
     */
    GPtrArray* solutions_pool = g_ptr_array_copy(
        root->localColPool, g_copy_scheduleset, &(problem.nb_jobs));

    problem.tree = std::make_unique<BranchBoundTree>(problem.root_pd, 0, 1);
    // call_branch_and_bound_explore(problem.tree);
    problem.tree->explore();

    /**
     * @brief Calculation of LB at the root node with column generation
     *
     */
    // if (problem.opt_sol->tw + problem.opt_sol->off != 0) {
    //     CCutil_start_timer(&(statistics->tot_lb_root));
    // compute_lower_bound(&problem, &(problem.root_pd));
    //     if (parms->pricing_solver < dp_solver) {
    //         solution_canonical_order(problem.opt_sol,
    //         root->local_intervals);
    //     }
    //     // represent_solution(root, problem.opt_sol);
    //     problem.rel_error =
    //         (double)(problem.global_upper_bound -
    //         problem.global_lower_bound) / (problem.global_lower_bound +
    //         0.00001);
    //     CCutil_stop_timer(&(statistics->tot_lb_root), 1);

    //     if (parms->pricing_solver == dp_bdd_solver) {
    //         int* take = get_take(root->solver);
    //         lp_node_data_free(root);
    //         root->localColPool = g_ptr_array_copy(
    //             solutions_pool, g_copy_scheduleset, &(problem.nb_jobs));
    //         build_rmp(root);
    //         freeSolver(root->solver);
    //         root->solver =
    //             newSolver(root->jobarray, root->nb_machines,
    //             root->ordered_jobs,
    //                       parms, problem.H_max, take,
    //                       problem.opt_sol->tw);

    //         CC_IFFREE(take, int);
    //         CCutil_start_timer(&(statistics->tot_lb_root));
    //         compute_lower_bound(&problem, root);
    //         CCutil_stop_timer(&(statistics->tot_lb_root), 1);
    //     }
    // }

    g_ptr_array_free(solutions_pool, TRUE);

    if (parms->mip_solver) {
        build_solve_mip(root);
    }

    if (parms->print) {
        problem.print_to_csv();
    }

    // CLEAN:
    // problem_free(&problem);
    return val;
}
