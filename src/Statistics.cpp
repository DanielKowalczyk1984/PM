
// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "Statistics.h"
#include <limits>
#include "Parms.h"         // for Parms
#include "orutils/util.h"  // for CCutil_init_timer, CCutil_start_resume_time, CCu...

Statistics::Statistics(const Parms& _parms)
    : global_upper_bound(std::numeric_limits<int>::max()),
      global_lower_bound(0),
      rel_error(std::numeric_limits<double>::max()),
      root_upper_bound(std::numeric_limits<int>::max()),
      root_lower_bound(0),
      root_rel_error(std::numeric_limits<double>::max()),
      nb_generated_col(0),
      nb_generated_col_root(0),
      first_size_graph(0),
      size_graph_after_reduced_cost_fixing(0),
      real_time_build_dd(0.0),
      real_time_total(getRealTime()),
      real_time_branch_and_bound(0.0),
      real_time_strong_branching(0.0),
      real_time_lb_root(0.0),
      real_time_lb(0.0),
      real_time_solve_lp(0.0),
      real_time_pricing(0.0),
      real_time_heuristic(0.0),
      mip_nb_vars(0),
      mip_nb_constr(0),
      mip_obj_bound(0.0),
      mip_obj_bound_lp(0.0),
      mip_rel_gap(0.0),
      mip_run_time(110.0),
      mip_status(0),
      mip_nb_iter_simplex(),
      mip_nb_nodes(0),
      mip_reduced_cost_fixing(),
      pname(_parms.pname) {
    CCutil_init_timer(&(tot_build_dd), "tot_build_dd");
    CCutil_init_timer(&(tot_cputime), "tot_cputime");
    CCutil_init_timer(&(tot_branch_and_bound), "tot_branch_and_bound");
    CCutil_init_timer(&(tot_strong_branching), "tot_strong_branching");
    CCutil_init_timer(&(tot_lb_root), "tot_lb_root");
    CCutil_init_timer(&(tot_lb), "tot_lb");
    CCutil_init_timer(&(tot_solve_lp), "tot_solve_lp");
    CCutil_init_timer(&(tot_pricing), "tot_pricing");
    CCutil_init_timer(&(tot_heuristic), "tot_heuristic");
    CCutil_init_timer(&(tot_reduce_cost_fixing), "tot_reduce_cost_fixing");
    CCutil_start_timer(&(tot_cputime));
}

void Statistics::start_resume_timer(TimerType _type) {
    switch (_type) {
        case build_dd_timer:
            CCutil_start_resume_time(&tot_build_dd);
            break;
        case cputime_timer:
            CCutil_start_resume_time(&tot_cputime);
            break;
        case bb_timer:
            CCutil_start_resume_time(&tot_branch_and_bound);
            break;
        case lb_root_timer:
            CCutil_start_resume_time(&tot_lb_root);
            break;
        case lb_timer:
            CCutil_start_resume_time(&tot_lb);
            break;
        case solve_lp_timer:
            CCutil_start_resume_time(&tot_solve_lp);
            break;
        case pricing_timer:
            CCutil_start_resume_time(&tot_pricing);
            break;
        case heuristic_timer:
            CCutil_start_resume_time(&tot_heuristic);
            break;
        case reduced_cost_fixing_timer:
            CCutil_start_resume_time(&tot_reduce_cost_fixing);
            break;
    }
}

void Statistics::suspend_timer(TimerType _type) {
    switch (_type) {
        case build_dd_timer:
            CCutil_suspend_timer(&tot_build_dd);
            break;
        case cputime_timer:
            CCutil_suspend_timer(&tot_cputime);
            break;
        case bb_timer:
            CCutil_suspend_timer(&tot_branch_and_bound);
            break;
        case lb_root_timer:
            CCutil_suspend_timer(&tot_lb_root);
            break;
        case lb_timer:
            CCutil_suspend_timer(&tot_lb);
            break;
        case solve_lp_timer:
            CCutil_suspend_timer(&tot_solve_lp);
            break;
        case pricing_timer:
            CCutil_suspend_timer(&tot_pricing);
            break;
        case heuristic_timer:
            CCutil_suspend_timer(&tot_heuristic);
            break;
        case reduced_cost_fixing_timer:
            CCutil_suspend_timer(&tot_reduce_cost_fixing);
            break;
    }
}

double Statistics::total_timer(TimerType _type) {
    switch (_type) {
        case build_dd_timer:
            return CCutil_total_timer(&tot_build_dd, 0);
        case cputime_timer:
            return CCutil_total_timer(&tot_cputime, 0);
        case bb_timer:
            return CCutil_total_timer(&tot_branch_and_bound, 0);
        case lb_root_timer:
            return CCutil_total_timer(&tot_lb_root, 0);
        case lb_timer:
            return CCutil_total_timer(&tot_lb, 0);
        case solve_lp_timer:
            return CCutil_total_timer(&tot_solve_lp, 0);
        case pricing_timer:
            return CCutil_total_timer(&tot_pricing, 0);
        case heuristic_timer:
            return CCutil_total_timer(&tot_heuristic, 0);
        case reduced_cost_fixing_timer:
            return CCutil_total_timer(&tot_reduce_cost_fixing, 0);
    }
    return 0.0;
}
