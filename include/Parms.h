// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef __PARMS_H__
#define __PARMS_H__

#include <fmt/format.h>
#include <array>                  // for array
#include <boost/timer/timer.hpp>  // for boost::timer::nanosecond_type
#include <cstddef>                // for size_t
#include <functional>             // for function
#include <nlohmann/json.hpp>      // for json
#include <string>                 // for allocator, string

enum BBNodeSelection {
    min_search_strategy = 0,
    no_branching = min_search_strategy,
    min_lb_strategy = 1,
    dfs_strategy = 2,
    max_strategy = 3,
};

enum PricingSolver {
    bdd_solver_simple = 0,
    bdd_solver_cycle = 1,
    bdd_solver_backward_simple = 2,
    bdd_solver_backward_cycle = 3,
    zdd_solver_simple = 4,
    zdd_solver_cycle = 5,
    zdd_solver_backward_simple = 6,
    zdd_solver_backward_cycle = 7,
    dp_solver = 8,
    ati_solver = 9,
    dp_bdd_solver = 10
};

enum StabTechniques {
    no_stab = 0,
    stab_wentgnes = 1,
    stab_dynamic = 2,
    stab_hybrid = 3,
    min_stab = stab_wentgnes,
};

enum BBExploreStrategy {
    min_bb_explore_strategy = 0,
    bb_dfs_strategy = min_bb_explore_strategy,
    bb_bfs_strategy = 1,
    bb_brfs_strategy = 2,
    bb_cbfs_strategy = 3,
};

enum Scoring_Parameter {
    min_scoring_parameter = 0,
    product_scoring_parameter = min_scoring_parameter,
    min_function_scoring_parameter = 1,
    weighted_sum_scoring_parameter = 2,
    weighted_product_scoring_parameter = 3,
    max_function_scoring_parameter = 4,
};

enum Scoring_Value {
    min_scoring_value = 0,
    lb_scoring_value = min_scoring_value,
    size_scoring_value = 1,
    nb_paths_scoring_value = 2,
};

NLOHMANN_JSON_SERIALIZE_ENUM(PricingSolver,
                             {{bdd_solver_simple, "BddForward"},
                              {bdd_solver_cycle, "BddForwardCycle"},
                              {bdd_solver_backward_simple, "BddBackward"},
                              {bdd_solver_backward_cycle, "BddBackwardCycle"},
                              {zdd_solver_simple, "ZddForward"},
                              {zdd_solver_cycle, "ZddForwardCycle"},
                              {zdd_solver_backward_simple, "ZddBackward"},
                              {zdd_solver_backward_cycle, "ZddBackwardCycle"},
                              {dp_solver, "Time-Indexed"},
                              {ati_solver, "Arc-Time-Indexed"},
                              {dp_bdd_solver, "Hybrid"}})

NLOHMANN_JSON_SERIALIZE_ENUM(StabTechniques,
                             {
                                 {no_stab, "NoStabilization"},
                                 {stab_wentgnes, "WentgnesStab"},
                                 {stab_dynamic, "DynamicStab"},
                                 {stab_hybrid, "HybridStab"},
                             })

NLOHMANN_JSON_SERIALIZE_ENUM(BBExploreStrategy,
                             {{min_bb_explore_strategy, "dfs"},
                              {bb_dfs_strategy, "dfs"},
                              {bb_bfs_strategy, "bfs"},
                              {bb_brfs_strategy, "brfs"},
                              {bb_cbfs_strategy, "cbfs"}})

NLOHMANN_JSON_SERIALIZE_ENUM(Scoring_Parameter,
                             {{min_scoring_parameter, "ProductScoring"},
                              {product_scoring_parameter, "ProductScoring"},
                              {min_function_scoring_parameter, "MinFunction"},
                              {max_function_scoring_parameter, "MaxFunction"},
                              {weighted_sum_scoring_parameter, " WeightedSum"},
                              {weighted_product_scoring_parameter,
                               "WeightedProduct"}})

struct Parms {
    template <typename T>
    class Parameter {
        std::string _name;
        T           _value;

       public:
        Parameter(const std::string name, const T& value)
            : _name(name),
              _value(value) {}

        inline const T&    value() const { return _value; }
        inline T&          value() { return _value; }
        const std::string& name() const { return _name; }
        std::string&       name() { return _name; }

        void set_value(const T& value) { _value = value; }
    };

    Parameter<BBExploreStrategy>                        bb_explore_strategy;
    Parameter<StabTechniques>                           stab_technique;
    Parameter<PricingSolver>                            pricing_solver;
    Parameter<Scoring_Parameter>                        scoring_parameter;
    Parameter<Scoring_Value>                            scoring_value;
    Parameter<int>                                      strong_branching;
    Parameter<int>                                      bb_node_limit;
    Parameter<int>                                      nb_iterations_rvnd;
    Parameter<boost::timer::nanosecond_type>            branching_cpu_limit;
    Parameter<bool>                                     use_cpu_time;
    Parameter<double>                                   alpha;
    Parameter<double>                                   branching_point;
    Parameter<bool>                                     use_heuristic;
    std::function<double(const std::array<double, 2>&)> scoring_function;
    Parameter<bool>                                     use_mip_solver;
    Parameter<bool>                                     refine_bdd;
    Parameter<bool>                                     enumerate;
    Parameter<bool>                                     pruning_test;
    Parameter<bool>                                     suboptimal_duals;
    Parameter<bool>                                     reduce_cost_fixing;
    Parameter<bool>                                     print_csv;

    /**
     * column generation
     */
    std::string jobfile;
    std::string pname;

    int    nb_jobs;
    size_t nb_machines;

    Parms();
    Parms(int argc, const char** argv);

    Parms(const Parms&) = default;
    Parms(Parms&&) = default;
    Parms& operator=(const Parms&) = default;
    Parms& operator=(Parms&&) = default;
    ~Parms() = default;

    void        parms_set_scoring_function(int scoring);
    void        parse_cmd(int argc, const char** argv);
    friend void from_json(const nlohmann::json& j, Parms& p);
    friend void to_json(nlohmann::json& j, const Parms& p);

   private:
    static constexpr double                mu = 5.0 / 6.0;
    static constexpr std::array<double, 2> beta = {1.5, 0.5};
    static constexpr double                TargetBrTimeValue = 0.45;
    static constexpr auto                  EPS = 1e-6;
};

template <typename T>
struct fmt::formatter<Parms::Parameter<T>> {
    char presentation = 'v';

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && (*it == 'v' || *it == 'n'))
            presentation = *it++;

        // Check if reached the end of the range:
        if (it != end && *it != '}')
            throw format_error("invalid format");

        return it;
    }

    template <typename FormatContext>
    auto format(const Parms::Parameter<T>& p, FormatContext& ctx)
        -> decltype(ctx.out()) {
        if (presentation == 'n') {
            return format_to(ctx.out(), "{}", p.name());
        }
        return format_to(ctx.out(), "{}", p.value());
    }
};

template <>
struct fmt::formatter<Parms> {
    char presentation = 'v';

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && (*it == 'v' || *it == 'n'))
            presentation = *it++;

        // Check if reached the end of the range:
        if (it != end && *it != '}')
            throw format_error("invalid format");

        return it;
    }

    template <typename FormatContext>
    auto format(const Parms& parms, FormatContext& ctx) -> decltype(ctx.out()) {
        auto aux_format =
            (presentation == 'n')
                ? "{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n}"
                : "{},{},{},{},{},{},{},{},{},{},{},{}";

        return format_to(
            ctx.out(), aux_format,
            parms.nb_iterations_rvnd, parms.stab_technique, parms.alpha,
            parms.pricing_solver, parms.strong_branching, parms.branching_point,
            parms.refine_bdd, parms.pruning_test, parms.suboptimal_duals,
            parms.scoring_parameter, parms.scoring_value, parms.use_cpu_time);
    }
};

#endif  // __PARMS_H__