////////////////////////////////////////////////////////////////
//                                                            //
//  main.c                                                    //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <defs.h>
#include <wct.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void usage(char *f) {
    fprintf(stderr, "Usage %s: [-see below-] adjlist_file NrMachines\n", f);
    fprintf(stderr, "   -d      turn on the debugging\n");
    fprintf(stderr,
            "   -f n    Number of feasible solutions that have to be "
            "constructed\n");
    fprintf(stderr,
            "   -s int  Node selection: 0 = none, 1= minimum lower "
            "bound(default), 2 = DFS\n");
    fprintf(stderr, "   -l dbl  Cpu time limit for branching\n");
    fprintf(stderr, "   -L dbl  Cpu time limit for scatter search\n");
    fprintf(stderr,
            "   -C int  Combine method scatter search: 0 = Pathrelinking, 1 = "
            "PM\n");
    fprintf(stderr,
            "   -r int  Scatter search use: 0 = no scatter "
            "search(default), 1 = scatter search\n");
    fprintf(stderr,
            "   -B int  Branch and Bound use: 0 = no branch and bound, 1 "
            "= use branch and bound(default)\n");
    fprintf(stderr,
            "   -S int  Stabilization technique: 0 = no "
            "stabilization(default), 1 = stabilization wentgnes, 2 = "
            "stabilization dynamic\n");
    fprintf(stderr,
            "   -z int  Pricing solver technique: 0 = BDD(default), 1 = "
            "ZDD, 2 = DP\n");
    fprintf(stderr,
            "   -c int  Construct heuristically solutions: 0 = "
            "yes(default), 1 = no\n");
    fprintf(stderr, "   -t int  Use ahv test: 0 = no(default), 1 = yes\n");
    fprintf(stderr, "   -p int  Print csv-files: 0 = no(default), 1 = yes\n");
    fprintf(stderr,
            "   -b int  Branching strategy: 0 = conflict(default), 1 = ahv\n");
    fprintf(stderr,
            "   -Z int  Use strong branching: 0 = use strong "
            "branching(default), 1 = no strong branching\n");
}

static int parseargs(int ac, char **av, wctparms *parms) {
    int c;
    int val = 0;
    int debug = dbg_lvl();

    while ((c = getopt(ac, av, "dr:f:s:l:L:C:B:z:S:c:D:t:p:b:Z:")) != EOF) {
        switch (c) {
            case 'd':
                ++(debug);
                set_dbg_lvl(debug);
                break;

            case 'r':
                val = wctparms_set_scatter_search(parms, atoi(optarg));
                CCcheck_val(val, "Failed cclasses_infile");
                break;

            case 'f':
                val = wctparms_set_nb_feas_sol(parms, atoi(optarg));
                CCcheck_val(val, "Failed number feasible solutions");
                break;

            case 's':
                val = wctparms_set_search_strategy(parms, atoi(optarg));
                CCcheck_val(val, "Failed set_branching_strategy");
                break;

            case 'l':
                val = wctparms_set_branching_cpu_limit(parms, atof(optarg));
                CCcheck_val(val, "Failed wctparms_set_branching_cpu_limit");
                break;

            case 'L':
                val =
                    wctparms_set_scatter_search_cpu_limit(parms, atof(optarg));
                CCcheck_val(val,
                            "Failed wctparms_set_scatter_search_cpu_limit");
                break;

            case 'C':
                val = wctparms_set_combine_method(parms, atoi(optarg));
                CCcheck_val(val, "Failed wctparms_set_combine_method");
                break;

            case 'B':
                val = wctparms_set_branchandbound(parms, atoi(optarg));
                CCcheck_val(val, "Failed wctparms_set_branchandbound");
                break;

            case 'S':
                val = wctparms_set_stab_technique(parms, atoi(optarg));
                CCcheck_val(val, "Failed in wctparms_set_stab_technique");
                break;

            case 'z':
                val = wctparms_set_solver(parms, atoi(optarg));
                CCcheck_val(val, "Failed in wctparms_set_solver");
                break;

            case 'c':
                val = wctparms_set_construct(parms, atoi(optarg));
                CCcheck_val(val, "Failed in construct sol");
                break;

            case 't':
                val = wctparms_set_test_ahv(parms, atoi(optarg));
                CCcheck_val(val, "Failed in use_test");
                break;

            case 'p':
                val = wctparms_set_print(parms, atoi(optarg));
                CCcheck_val(val, "Failed in print");
                break;

            case 'b':
                val = wctparms_set_branching_strategy(parms, atoi(optarg));
                CCcheck_val(val, "Failed in set branching strategy");
                break;

            case 'Z':
                val = wctparms_set_strong_branching(parms, atoi(optarg));
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
        val = wctparms_set_file(parms, av[optind++]);
        CCcheck_val(val, "Failed in wctparms_set_file");

        if (ac <= optind) {
            val = 1;
            goto CLEAN;
        }

        val = wctparms_set_nmachines(parms, atoi(av[optind++]));
        CCcheck_val(val, "Failed in wctparms_set_nmachines");
    }

CLEAN:

    if (val) {
        usage(av[0]);
    }

    return val;
}

int main(int ac, char **av) {
    int           val = 0;
    double start_time;
    wctproblem    problem;
    wctdata *     pd;
    wctparms *    parms;

    val = program_header(ac, av);
    CCcheck_val_2(val, "Failed in program_header")
    wctproblem_init(&problem);
    CCutil_start_timer(&(problem.tot_cputime));
    start_time = CCutil_zeit();
    problem.real_time = getRealTime();
    problem.nwctdata = 1;
    pd = &(problem.root_pd);
    wctdata_init(pd);
    pd->id = 0;
    parms = &(problem.parms);
    val = parseargs(ac, av, parms);
    CCcheck_val_2(val, "Failed in parseargs");

    if (dbg_lvl() > 1) {
        printf("Debugging turned on\n");
    }

    fflush(stdout);
    /** Reading and preprocessing the data */
    val = read_problem(&problem);
    CCcheck_val_2(val, "read_adjlist failed");
    val = preprocess_data(&problem);
    CCcheck_val_2(val, "Failed at preprocess_data");

    /** Finding heuristic solutions to the problem */
    heuristic_rpup(&problem);
    printf("Reading and preprocessing of the data took %f seconds\n",
           CCutil_zeit() - start_time);

    /** Branch-and-Price Algorithm */

CLEAN:
    wctproblem_free(&problem);
    return val;
}
