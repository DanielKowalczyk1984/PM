#ifndef __USAGE_H__
#define __USAGE_H__

#include <string>

static const char USAGE[] = R"(PM.
  Usage: 
    bin/PM FILE NB --json <json_file>
    bin/PM --version
    bin/PM (-h | --help)

  Arguments:
    FILE      Path to the instance file
    NB        Number of machines

  options:
    -h --help                     Show this screen.
    --version                     Show version.
    -j --json                     Use json settings file.
)";
    // bin/PM FILE NB [-s <sn> -S <kn> -pRZHMdPD -n <nl> -b <br> -a <ln> -l <x> -f <y> -c <x> --alpha <mn> --branching_point <brp> --refinement --enumerate --scoring_value <sv>]
    // -S --stab_method=<kn>         Stabilization technique: 0 = no stabilization, 1 = stabilization wentgnes, 2 = stabilization dynamic[default: 1].
    // -l --cpu_limit=<x>            Cpu time limit for branch and bound method[default: 7200].
    // -f --nb_rvnb_it=<y>           Number of iterations in RVND[default: 5000].
    // --alpha=<mn>                  Stabilization factor[default: 0.8].
    // --scoring_value=<sv>          Scoring value[default:0].  
    // --branching_point=<brp>       Branching point[default: 0.5].
    // -p --print_csv                Print csv-files.
    // -r --refinement               Refine decision diagram.
    // -e --enumerate                Enumerate elementary paths.
    // -R --no_rc_fixing             Don't apply reduce cost fixing.
    // -H --no_heuristic             Don't apply heuristic.
    // -c --strong_branching=<sb>    Don't apply strong branching[default: 8].
    // -M --use_mip_solver           Use MIP solver.
    // -a --pricing_solver=<ln>      Set pricing solver: 0 = bdd backward cycle, 1 = bdd forward simple, 2 = bdd forward cycle, 3 = bdd backward simple, 4 = bdd backward cycle, 5 = zdd forward simple, 6 = zdd forward cycle, 7 = zdd backward simple, 8 = zdd backward cycle, 9 = TI solver, 10 = arc-TI solver, 11 = hybrid model TI and bdd backward cycle[default: 4].
    // -b --branching_strategy=<br>  Set branch-and-bound exploration strategy: 0 = DFS, 1 = BFS, 2 = BrFS, 3 = CBFS[default: 0].
    // -n --node_limit=<nl>          Set a limit on the number of nodes that can be explored.[default: 0]. Default meaning that all nodes will be explored.
    // -P --pruning_test             Use pruning test in branch and bound tree.
    // -D --suboptimal_duals         Use sub_optimal dual variables.
    // -d --debug                    Turn on the debugging.
    // -s --scoring_function=<sn>    Set scoring function branching[default: 0]
#endif  // __USAGE_H__