#include <wctparms.h>
#include <regex>
#include <string>

#include <docopt/docopt.h>

static const std::string USAGE =
    R"(PM.

Usage:
  bin/PM [-S <kn> -pmBRZH -b <br> -a <ln> -l <x> -f <y> -d --no_strong_branching --alpha <mn>] FILE NB
  bin/PM (-h | --help)
  bin/PM --version

Arguments:
  FILE  Path to the instance file
  NB    Number of machines

Options:
  -h --help                     Show this screen.
  --version                     Show version.
  -d --debug                    Turn on the debugging.
  -S --stab_method=<kn>         Stabilization technique: 0 = no stabilization, 1 = stabilization wentgnes, 2 = stabilization dynamic[default: 1].
  -a --pricing_solver=<ln>      Set pricing solver: 0 = bdd backward cycle, 1 = bdd forward simple, 2 = bdd forward cycle,
                                                    3 = bdd backward simple, 4 = bdd backward cycle, 5 = zdd forward simple,
                                                    6 = zdd forward cycle, 7 = zdd backward simple, 8 = zdd backward cycle,
                                                    9 = TI solver, 10 = arc-TI solver, 11 = hybrid model TI and bdd backward cycle[default: 4].
  -l --cpu_limit=<x>            Cpu time limit for branch and bound method[default: 7200].
  -f --nb_rvnb_it=<y>           Number of iterations in RVND[default: 4].
  --alpha=<mn>                  Stabilization factor[default: 0.8].
  -p --print_csv                Print csv-files.
  -m --mip_solver               Use mip solver to solve the original formulation.
  -R --no_rc_fixing             Don't apply reduce cost fixing.
  -H --no_heuristic             Don't apply heuristic.
  -B --no_branch_and_bound      Don't apply branch-and-bound.
  --no_strong_branching         Don't apply strong branching.
  -b --branching_strategy=<br>  Set branch-and-bound exploration strategy: 0 = DFS, 1 = BFS, 2 = BrFS, 3 = CBFS[default: 0].
)";

static std::string find_match(const char* _instance_file) {
    std::regex  regexp{"^.*(wt[0-9]*_[0-9]*).*$"};
    std::smatch match{};
    std::string s{_instance_file};
    std::regex_search(s, match, regexp);

    if (match.size() != 2) {
        return std::string("unknown_problem");
    } else {
        return match[1];
    }
}

extern "C" int parse_cmd(int argc, const char** argv, Parms* parms) {
    int val = 0;

    auto args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, "PM 0.1");

    /** Set CPU limit for branching */
    val = parms_set_branching_cpu_limit(
        parms, static_cast<double>(args["--cpu_limit"].asLong()));
    /** Set number of iterations in rvnd */
    val = parms_set_nb_iterations_rvnd(
        parms, static_cast<int>(args["--nb_rvnb_it"].asLong()));
    /** Print statistics to csv files */
    val = parms_set_print(parms, args["--print_csv"].asBool());
    /** Use MIP solver or not */
    val = parms_set_mip_solver(parms, args["--mip_solver"].asBool());
    /** Use reduced cost fixing */
    val = parms_set_reduce_cost(parms, !(args["--no_rc_fixing"].asBool()));
    /** Use heuristic or not */
    val = parms_set_use_heuristic(parms, !(args["--no_heuristic"].asBool()));
    /** Set the pricing solver */
    val = parms_set_pricing_solver(
        parms, static_cast<int>(args["--pricing_solver"].asLong()));
    /** Set the stabilization method */
    val = parms_set_stab_technique(
        parms, static_cast<int>(args["--stab_method"].asLong()));
    val =
        parms_set_branchandbound(parms, args["--no_branch_and_bound"].asBool());
    /** Determine the name of the instance */
    auto* file_name = args["FILE"].asString().c_str();
    val = parms_set_file(parms, file_name);
    val = parms_set_pname(parms, find_match(file_name).c_str());
    /** Set the number of machines */
    val = parms_set_nb_machines(parms, static_cast<int>(args["NB"].asLong()));

    val = parms_set_strong_branching(parms,
                                     !(args["--no_strong_branching"].asBool()));
    val = parms_set_alpha(parms, std::stod(args["--alpha"].asString()));
    val = parms_set_bb_explore_strategy(
        parms, static_cast<int>(args["--branching_strategy"].asLong()));
    return val;
}
