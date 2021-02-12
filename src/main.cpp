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

    problem.solve();
    // CLEAN:
    // problem_free(&problem);
    return val;
}
