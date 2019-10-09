
////////////////////////////////////////////////////////////////
//                                                            //
//  IP_models.c                                                    //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 08/10/2019.                            //
//  Copyright (c) 2019 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <defs.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <wct.h>
#include <wctparms.h>
static void usage(char* f) {
    fprintf(stderr, "Usage %s: [-see below-] adjlist_file NrMachines\n", f);
    fprintf(stderr, "   -d      turn on the debugging\n");
    fprintf(stderr, "   -f n    Number of iterations in RVND\n");
    fprintf(stderr,
            "   -s int  Node selection: 0 = none, 1= minimum lower "
            "bound(default), 2 = DFS\n");
    fprintf(stderr, "   -l dbl  Cpu time limit for branching\n");
    fprintf(stderr, "   -L dbl  Cpu time limit for heuristic construction\n");
    fprintf(stderr,
            "   -B int  Branch and Bound use: 0 = no branch and bound, 1 "
            "= use branch and bound(default)\n");
    fprintf(
        stderr,
        "   -S int  Stabilization technique: 0 = no stabilization(default), "
        "1 = stabilization wentgnes, 2 = stabilization dynamic\n");
    fprintf(stderr, "   -p int  Print csv-files: 0 = no(default), 1 = yes\n");
    fprintf(stderr,
            "   -b int  Branching strategy: 0 = conflict(default), 1 = ahv\n");
    fprintf(stderr,
            "   -Z int  Use strong branching: 0 = use strong "
            "branching(default), 1 = no strong branching\n");
    fprintf(stderr,
            "   -a int  Set solver: 0 = bdd solver(default), 1 = dp solver\n");
}

static char* find_match(const char* _instance_file) {
    regex_t    regex;
    regmatch_t matches[2];
    char* val = (char*) NULL;
    int regex_val = regcomp(&regex, "^.*(wt[0-9]*_[0-9]*).*$", REG_EXTENDED);
    if (regex_val) {
        fprintf(stderr, "Could not compile regex\n");
        val = strdup("regex_error");
    } else {
        regex_val = regexec(&regex, _instance_file, 2, (regmatch_t*)&matches, 0);
        if (!regex_val) {
            val = g_strndup(_instance_file + matches[1].rm_so,
                                matches[1].rm_eo - matches[1].rm_so);
        } else {
            printf("No match\n");
            val = strdup("unknown_name");
        }
    }
    regfree(&regex);

    return val;
}

static int parseargs(int ac, char** av, Parms* parms) {
    int    c;
    double f;
    int    val = 0;
    char*  ptr;
    int    debug = dbg_lvl();

    while ((c = getopt(ac, av, "df:s:l:L:B:S:D:p:b:Z:a:")) != EOF) {
        switch (c) {
            case 'd':
                ++(debug);
                set_dbg_lvl(debug);
                break;

            case 'f':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_nb_iterations_rvnd(parms, c);
                CCcheck_val(val, "Failed number feasible solutions");
                break;

            case 's':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_search_strategy(parms, c);
                CCcheck_val(val, "Failed set_branching_strategy");
                break;

            case 'l':
                f = strtod(optarg, &ptr);
                val = parms_set_branching_cpu_limit(parms, f);
                CCcheck_val(val, "Failed parms_set_branching_cpu_limit");
                break;

            case 'L':
                f = strtod(optarg, &ptr);
                val = parms_set_scatter_search_cpu_limit(parms, f);
                CCcheck_val(val, "Failed parms_set_scatter_search_cpu_limit");
                break;

            case 'B':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_branchandbound(parms, c);
                CCcheck_val(val, "Failed parms_set_branchandbound");
                break;

            case 'S':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_stab_technique(parms, c);
                CCcheck_val(val, "Failed in parms_set_stab_technique");
                break;

            case 'p':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_print(parms, c);
                CCcheck_val(val, "Failed in print");
                break;

            case 'b':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_branching_strategy(parms, c);
                CCcheck_val(val, "Failed in set branching strategy");
                break;
            case 'a':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_pricing_solver(parms, c);
                CCcheck_val(val, "Failed in set alpha");
                break;

            case 'Z':
                c = strtol(optarg, &ptr, 10);
                val = parms_set_strong_branching(parms, c);
                CCcheck_val(val, "Failed in set strong branching");
                break;

            default:
                usage(av[0]);
                val = 1;
                goto CLEAN;
        }
    }

    if (ac <= optind) {
        val = 1;
        goto CLEAN;
    } else {
        char* instance_file = av[optind++];

        val = parms_set_file(parms, instance_file);
        CCcheck_val(val, "Failed in parms_set_file");
        val = parms_set_pname(parms, find_match(instance_file));
        CCcheck_val(val, "Failed in parms_set_pname");

        if (ac <= optind) {
            val = 1;
            goto CLEAN;
        }
        c = strtol(av[optind++], &ptr, 10);
        val = parms_set_nb_machines(parms, c);
        CCcheck_val(val, "Failed in parms_set_nb_machines");
    }

CLEAN:

    if (val) {
        usage(av[0]);
    }

    return val;
}

// int main(int ac, char** av) {
//     int       val = 0;
//     double    start_time;
//     Problem   problem;
//     NodeData* root = &(problem.root_pd);
//     Parms*    parms = &(problem.parms);
//     val = program_header(ac, av);
//     CCcheck_val_2(val, "Failed in program_header");
//     problem_init(&problem);
//     val = parseargs(ac, av, &(problem.parms));
//     CCcheck_val_2(val, "Failed in parseargs");
//     if (dbg_lvl() > 1) {
//         printf("Debugging turned on\n");
//     }

//     /**
//      * Reading and preprocessing the data
//      */
//     start_time = CCutil_zeit();
//     val = read_problem(&problem);
//     CCcheck_val_2(val, "read_adjlist failed");
//     val = preprocess_data(&problem);
//     CCcheck_val_2(val, "Failed at preprocess_data");
//     printf("Reading and preprocessing of the data took %f seconds\n",
//            CCutil_zeit() - start_time);

//     /**
//      * Finding heuristic solutions to the problem
//      */
//     heuristic(&problem);

//     /**
//      * Build DD at the root node
//      */
//     if (parms->pricing_solver < dp_solver) {
//         CCutil_start_timer(&(problem.tot_build_dd));
//         root->solver = newSolver(root->jobarray, root->nb_machines,
//                                  root->ordered_jobs, &(problem.parms));
//         CCutil_stop_timer(&(problem.tot_build_dd), 0);
//         print_size_to_csv(&problem, root);
//     } else {
//         root->solver =
//             newSolverDp(root->jobarray, root->nb_machines, root->H_max, parms);
//     }
//     problem.first_size_graph = get_nb_edges(root->solver);

//     /**
//      * Calculation of LB at the root node with column generation
//      */
//     if (problem.opt_sol->tw + problem.opt_sol->off != 0) {
//         build_rmp(&(problem.root_pd), 0);
//         CCutil_start_timer(&(problem.tot_lb_root));
//         compute_lower_bound(&problem, &(problem.root_pd));
//         print_dot_file(root->solver, NULL);
//         problem.rel_error =
//             (double)(problem.global_upper_bound - problem.global_lower_bound) /
//             (problem.global_lower_bound + 0.00001);
//         CCutil_stop_timer(&(problem.tot_lb_root), 1);
//         // represent_solution(root, problem.opt_sol);
//         // construct_lp_sol_from_rmp(root);
//         // build_solve_mip(root);
//         // build_solve_mip(root);
//         // represent_solution(root, problem.opt_sol);
//         // compute_lower_bound(&problem, &(problem.root_pd));
//         // problem.size_graph_after_reduced_cost_fixing =
//         // get_size_graph(root->solver);
//         // evaluate_nodes(root);
//         // lpwctdata_free(&(problem.root_pd));
//         // problem.root_pd.localColPool =
//         // g_ptr_array_new_with_free_func(g_scheduleset_free);
//         // heuristic_rpup(&problem);
//         // build_rmp(&(problem.root_pd), 0);
//         // CCutil_start_timer(&(problem.tot_lb_root));
//         // compute_lower_bound(&problem, &(problem.root_pd));
//         // CCutil_stop_timer(&(problem.tot_lb_root), 1);
//     }
//     print_to_csv(&problem);

// CLEAN:
//     problem_free(&problem);
//     return val;
// }
