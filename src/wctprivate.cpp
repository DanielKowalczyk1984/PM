#include "wctprivate.h"
#include <fmt/core.h>
#include <algorithm>
#include <memory>
#include <type_traits>
#include "PricerSolverArcTimeDP.hpp"
#include "PricerSolverBddBackward.hpp"
#include "PricerSolverBddForward.hpp"
#include "PricerSolverSimpleDP.hpp"
#include "PricerSolverZddBackward.hpp"
#include "PricerSolverZddForward.hpp"
#include "PricingStabilization.hpp"
// #include "interval.h"

Problem::~Problem() { /*free the parameters*/
    g_ptr_array_free(g_job_array, TRUE);
    g_ptr_array_free(ColPool, TRUE);
    g_ptr_array_free(intervals, TRUE);
    solution_free(&(opt_sol));
}

Problem::Problem(int argc, const char** argv) {
    double start_time = 0.0;
    problem_init();

    int val = program_header(argc, argv);
    val = parms.parse_cmd(argc, argv);

    if (dbg_lvl() > 1) {
        fmt::print("Debugging turned on\n");
    }

    /**
     * @brief Reading and preprocessing the data
     *
     */
    start_time = CCutil_zeit();
    problem_read();
    preprocess_data();
    // CCcheck_val_2(val, "Failed at preprocess_data");
    fmt::print("Reading and preprocessing of the data took %f seconds\n",
               CCutil_zeit() - start_time);
    /**
     *@brief Finding heuristic solutions to the problem or start without
     *feasible solutions
     *
     */
    if (parms.use_heuristic) {
        heuristic();
    } else {
        opt_sol = solution_alloc(intervals->len, nb_machines, nb_jobs, off);
        Solution* sol = opt_sol;
        // CCcheck_NULL_2(sol, "Failed to allocate memory");
        val = construct_edd(sol);
        // CCcheck_val_2(val, "Failed construct edd");
        fmt::print("Solution Constructed with EDD heuristic:\n");
        solution_print(sol);
        solution_canonical_order(sol, intervals);
        fmt::print("Solution in canonical order: \n");
        solution_print(sol);
    }
    /**
     * @brief Build DD at the root node
     *
     */
    CCutil_start_timer(&(stat.tot_build_dd));
    switch (parms.pricing_solver) {
        case bdd_solver_simple:
            root_pd->solver = std::make_unique<PricerSolverBddSimple>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), H_max, nullptr, opt_sol->tw);
            break;
        case bdd_solver_cycle:
            root_pd->solver = std::make_unique<PricerSolverBddCycle>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), H_max, nullptr, opt_sol->tw);
            break;
        case zdd_solver_cycle:
            root_pd->solver = std::make_unique<PricerSolverZddCycle>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), global_upper_bound);
            break;
        case zdd_solver_simple:
            root_pd->solver = std::make_unique<PricerSolverSimple>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), opt_sol->tw);
            break;
        case bdd_solver_backward_simple:
            root_pd->solver = std::make_unique<PricerSolverBddBackwardSimple>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), H_max, nullptr, global_upper_bound);
            break;
        case bdd_solver_backward_cycle:
            root_pd->solver = std::make_unique<PricerSolverBddBackwardCycle>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), H_max, nullptr, global_upper_bound);
            ;
            break;
        case zdd_solver_backward_simple:
            root_pd->solver = std::make_unique<PricerSolverZddBackwardSimple>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), opt_sol->tw);
            break;
        case zdd_solver_backward_cycle:
            root_pd->solver = std::make_unique<PricerSolverZddBackwardCycle>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), opt_sol->tw);
        case dp_solver:
            root_pd->solver = std::make_unique<PricerSolverSimpleDp>(
                g_job_array, nb_machines, H_max, parms.pname.c_str(),
                opt_sol->tw);
            break;
        case ati_solver:
            root_pd->solver = std::make_unique<PricerSolverArcTimeDp>(
                g_job_array, nb_machines, H_max, parms.pname.c_str(),
                opt_sol->tw);
            break;
        case dp_bdd_solver:
            root_pd->solver = std::make_unique<PricerSolverSimpleDp>(
                g_job_array, nb_machines, H_max, parms.pname.c_str(),
                opt_sol->tw);
            break;
        default:
            root_pd->solver = std::make_unique<PricerSolverBddBackwardCycle>(
                g_job_array, nb_machines, root_pd->ordered_jobs,
                parms.pname.c_str(), H_max, nullptr, global_upper_bound);
            break;
    }
    CCutil_stop_timer(&(stat.tot_build_dd), 0);
    stat.first_size_graph = root_pd->solver->get_nb_vertices();

    /**
     * @brief Initial stabilization method
     *
     */
    auto* tmp_solver = root_pd->solver.get();
    switch (parms.stab_technique) {
        case stab_wentgnes:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationStat>(tmp_solver);
            break;
        case stab_dynamic:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationDynamic>(tmp_solver);
            break;
        case stab_hybrid:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationHybrid>(tmp_solver);
            break;
        case no_stab:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationBase>(tmp_solver);
            break;
        default:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationStat>(tmp_solver);
            break;
    }
    root_pd->stat = &(stat);
    root_pd->opt_sol = opt_sol;

    /**
     * @brief Initialization of the B&B tree
     *
     */
    tree = std::make_unique<BranchBoundTree>(std::move(root_pd), 0, 1);
    tree->explore();
}

static void g_problem_summary_init(gpointer data, gpointer user_data) {
    Problem* prob = static_cast<Problem*>(user_data);
    Job*     j = static_cast<Job*>(data);

    prob->p_sum += j->processing_time;
    prob->pmax = std::max(prob->pmax, j->processing_time);
    prob->pmin = std::min(prob->pmin, j->processing_time);
    prob->dmax = std::max(prob->dmax, j->due_time);
    prob->dmin = std::min(prob->dmin, j->due_time);
}

void Problem::calculate_Hmax() {
    int       temp = 0;
    double    temp_dbl = 0.0;
    NodeData* pd = root_pd.get();

    temp = p_sum - pmax;
    temp_dbl = static_cast<double>(temp);
    temp_dbl = floor(temp_dbl / nb_machines);
    H_max = pd->H_max = static_cast<int>(temp_dbl) + pmax;
    H_min = pd->H_min = static_cast<int>(ceil(temp_dbl / nb_machines)) - pmax;

    GPtrArray* duration = g_ptr_array_copy(g_job_array, NULL, NULL);
    g_ptr_array_set_free_func(duration, NULL);
    g_ptr_array_sort(duration, g_compare_duration);

    int    m = 0;
    int    i = nb_jobs - 1;
    double tmp = p_sum;
    H_min = p_sum;
    do {
        Job* job = static_cast<Job*>(g_ptr_array_index(duration, i));
        tmp -= job->processing_time;
        m++;
        i--;

    } while (m < nb_machines - 1);

    H_min = pd->H_min = static_cast<int>(ceil(tmp / nb_machines));
    g_ptr_array_free(duration, TRUE);
    fmt::print(
        "H_max = {}, H_min = {},  pmax = {}, pmin = {}, p_sum = {}, off = "
        "{}\n",
        H_max, H_min, pmax, pmin, p_sum, off);
}

void Problem::create_ordered_jobs_array(GPtrArray* a, GPtrArray* b) {
    // interval* tmp_interval = (interval*)NULL;
    // Job*               tmp_j = nullptr;
    // job_interval_pair* tmp_pair = static_cast<job_interval_pair*)NULL;
    for (unsigned i = 0; i < a->len; ++i) {
        auto* tmp_interval = static_cast<interval*>(g_ptr_array_index(a, i));
        GPtrArray* jobarray = tmp_interval->sigma;
        for (unsigned j = 0; j < jobarray->len; ++j) {
            auto* tmp_j = static_cast<Job*>(g_ptr_array_index(jobarray, j));
            if (tmp_j->processing_time <= tmp_interval->b) {
                job_interval_pair* tmp_pair =
                    CC_SAFE_MALLOC(1, job_interval_pair);
                tmp_pair->j = tmp_j;
                tmp_pair->I = tmp_interval;
                g_ptr_array_add(b, tmp_pair);
            }
        }
    }

    fmt::print("There are {} layers\n", b->len);
}

int Problem::find_division() {
    int            val = 0;
    int            counter = 0;
    int            prev = 0;
    GPtrArray*     tmp_array = g_ptr_array_new_with_free_func(g_interval_free);
    Job*           tmp_j = NULL;
    Job *          j1 = NULL, *j2 = NULL;
    interval*      tmp_interval = NULL;
    interval_pair* pair = NULL;
    interval_pair  tmp_pair;

    /** Find initial partition */
    prev = 0;
    for (int i = 0; i < nb_jobs && prev < H_max; ++i) {
        tmp_j = static_cast<Job*>(g_ptr_array_index(g_job_array, i));
        int tmp = CC_MIN(H_max, tmp_j->due_time);
        if (prev < tmp) {
            tmp_interval = interval_alloc(prev, tmp, -1, g_job_array, nb_jobs);
            g_ptr_array_add(tmp_array, tmp_interval);
            CCcheck_NULL_2(tmp_interval, "Failed to allocate memory");
            prev = tmp_j->due_time;
        }
    }

    if (prev < H_max) {
        tmp_interval = interval_alloc(prev, H_max, -1, g_job_array, nb_jobs);
        g_ptr_array_add(tmp_array, tmp_interval);
    }

    /** calculate the new intervals */
    for (unsigned i = 0; i < tmp_array->len; ++i) {
        GList* pairs = (GList*)NULL;
        tmp_interval = static_cast<interval*>(g_ptr_array_index(tmp_array, i));
        for (size_t j = 0; j < tmp_interval->sigma->len - 1; j++) {
            for (size_t k = j + 1; k < tmp_interval->sigma->len; k++) {
                j1 = static_cast<Job*>(
                    g_ptr_array_index(tmp_interval->sigma, j));
                j2 = static_cast<Job*>(
                    g_ptr_array_index(tmp_interval->sigma, k));
                tmp_pair =
                    (interval_pair){j1, j2, tmp_interval->a, tmp_interval->b};
                if (!check_interval(&tmp_pair, i, tmp_array) &&
                    !(j1->due_time >= tmp_interval->b) &&
                    !(j2->due_time >= tmp_interval->b)) {
                    pair = CC_SAFE_MALLOC(1, interval_pair);
                    *pair =
                        (interval_pair){j1, j2, tmp_pair.left, tmp_pair.right};
                    pairs = g_list_append(pairs, pair);
                }
            }
        }

        if (pairs) {
            GPtrArray* slots = array_time_slots(tmp_interval, pairs);
            for (unsigned j = 1; j < slots->len; ++j) {
                g_ptr_array_add(
                    intervals, interval_alloc(*((int*)slots->pdata[j - 1]),
                                              *((int*)slots->pdata[j]), counter,
                                              g_job_array, nb_jobs));
                counter++;
            }
            g_ptr_array_free(slots, TRUE);
        } else {
            g_ptr_array_add(intervals,
                            interval_alloc(tmp_interval->a, tmp_interval->b,
                                           counter, g_job_array, nb_jobs));
            counter++;
        }
    }

CLEAN:
    g_ptr_array_free(tmp_array, TRUE);
    return val;
}

int Problem::check_interval(interval_pair* pair,
                            int            k,
                            GPtrArray*     interval_array) {
    auto* I = static_cast<interval*>(g_ptr_array_index(interval_array, k));
    auto* j = pair->b;
    return (I->a + j->processing_time >= I->b ||
            calculate_T(pair, k, interval_array) <= I->a);
}

int Problem::calculate_T(interval_pair* pair,
                         int            k,
                         GPtrArray*     interval_array) {
    auto* I = static_cast<interval*>(g_ptr_array_index(interval_array, k));
    auto* i = pair->a;
    auto* j = pair->b;
    pair->left = I->a;
    pair->right = I->a + j->processing_time;

    if (pair->left > I->b - i->processing_time) {
        return pair->left;
    } else {
        if (value_diff_Fij(pair->left, i, j) <= 0) {
            return pair->left;
        }

        // for (int t = k - 1; t >= 0; t--) {
        //     tmp = (interval*)g_ptr_array_index(interval_array, t);
        //     pair->left = tmp->a + j->processing_time - i->processing_time;

        //     if (value_diff_Fij(pair->left, i, j) <= 0 && pair->left >= tmp->a
        //     &&
        //         pair->left <= tmp->b - i->processing_time) {
        //         break;
        //     }
        // }

        pair->left =
            i->due_time +
            (int)ceil((double)(j->weight * i->processing_time) / i->weight) -
            i->processing_time;
        return pair->left;
    }
}

GPtrArray* Problem::array_time_slots(interval* I, GList* pairs) {
    GPtrArray* array = g_ptr_array_new_with_free_func(free);
    // interval_pair* tmp = (interval_pair*)NULL;
    interval_pair* min_data = nullptr;
    GList*         min = (GList*)NULL;
    int*           tmp_int = NULL;

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->a;
    g_ptr_array_add(array, tmp_int);

    while (pairs) {
        min = pairs;
        min_data = static_cast<interval_pair*>(min->data);
        for (GList* i = min->next; i; i = g_list_next(i)) {
            auto* tmp = static_cast<interval_pair*>(i->data);
            if (tmp->right < min_data->right) {
                min = i;
                min_data = (interval_pair*)i->data;
            }
        }

        tmp_int = CC_SAFE_MALLOC(1, int);
        *tmp_int = min_data->right;
        g_ptr_array_add(array, tmp_int);
        pairs = g_list_remove_link(pairs, min);
        g_list_free_full(min, interval_pair_free);

        GList* i = pairs;
        while (i) {
            auto* tmp = static_cast<interval_pair*>(i->data);
            if ((*tmp_int >= tmp->left && *tmp_int <= tmp->right) ||
                (*tmp_int + tmp->b->processing_time >= I->b)) {
                GList* remove = i;
                i = g_list_next(i);
                pairs = g_list_remove_link(pairs, remove);
                g_list_free_full(remove, interval_pair_free);
            } else {
                tmp->right = *tmp_int + tmp->b->processing_time;
                i = g_list_next(i);
            }
        }
    }

    tmp_int = CC_SAFE_MALLOC(1, int);
    *tmp_int = I->b;
    g_ptr_array_add(array, tmp_int);

    return array;
}

NodeData::~NodeData() {
    temporary_data_free();
    g_ptr_array_free(ordered_jobs, TRUE);
    g_ptr_array_free(best_schedule, TRUE);
}

int Problem::preprocess_data() {
    int val = 0;
    int i = 0;
    // NodeData* root = root_pd.get();

    /** Calculate the statistics of the instance */
    g_ptr_array_foreach(g_job_array, g_problem_summary_init, this);

    /** Calculate H_max */
    calculate_Hmax();

    /** order the jobarray of problem following edd rule */
    g_ptr_array_sort_with_data(g_job_array, g_job_compare_edd, NULL);

    g_ptr_array_foreach(g_job_array, g_set_jobarray_job, &i);
    root_pd->jobarray = g_job_array;
    root_pd->off = off;

    /** Find the intervals of the instance at the root node */
    find_division();

    /** Create all node of the ZDD */
    create_ordered_jobs_array(intervals, root_pd->ordered_jobs);

    /** Determine the position of each job in the interval */
    // determine_jobs_order_interval(problem);

    return val;
}

void Problem::solve() {
    tree->explore();
    if (parms.print) {
        print_to_csv();
    }
}
