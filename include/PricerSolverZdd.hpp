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

#ifndef PRICER_SOLVER_ZDD_HPP
#define PRICER_SOLVER_ZDD_HPP
#include <cstddef>                        // for size_t
#include <memory>                         // for unique_ptr, shared_ptr
#include <utility>                        // for pair
#include <vector>                         // for vector
#include "Instance.h"                     // for Instance
#include "MipGraph.hpp"                   // for MipGraph
#include "ModernDD/NodeBddStructure.hpp"  // for DdStructure
#include "PricerSolverBase.hpp"  // for PricerSolverBase, PricerSolverBase::...
#include "PricingSolution.hpp"   // for PricingSolution
#include "ZddNode.hpp"           // for NodeZdd
struct Column;                   // lines 19-19
struct Interval;                 // lines 16-16
struct Job;                      // lines 17-17
struct NodeData;                 // lines 18-18
class PricerSolverZdd : public PricerSolverBase {
   public:
    std::unique_ptr<DdStructure<NodeZdd<double>>> decision_diagram;

    size_t size_graph{};
    size_t nb_removed_edges{};
    size_t nb_removed_nodes{};

    // GPtrArray*                              ordered_jobs;
    std::vector<std::pair<Job*, Interval*>> ordered_jobs_new;

    MipGraph            mip_graph;
    std::vector<double> lp_x;
    std::vector<double> solution_x;

    explicit PricerSolverZdd(const Instance& instance);
    PricerSolverZdd(PricerSolverZdd&&) = default;
    PricerSolverZdd& operator=(PricerSolverZdd&&) = default;
    PricerSolverZdd& operator=(const PricerSolverZdd&) = delete;

    PricerSolverZdd(const PricerSolverZdd& src)
        : PricerSolverBase(src),
          size_graph(src.size_graph),
          nb_removed_edges(src.nb_removed_edges),
          nb_removed_nodes(src.nb_removed_nodes),
          ordered_jobs_new(src.ordered_jobs_new),
          mip_graph(src.mip_graph) {}

    ~PricerSolverZdd() override = default;
    [[nodiscard]] std::unique_ptr<PricerSolverBase> clone() const override {
        return nullptr;
    };

    void init_table();

    void remove_layers();
    void remove_edges();

    void construct_mipgraph();
    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const double*                               lambda,
        const std::vector<std::shared_ptr<Column>>& columns) override;
    void construct_lp_sol_from_rmp(
        const std::span<const double>&              lambda,
        const std::vector<std::shared_ptr<Column>>& columns) override;
    bool    check_column(Column const* set) override;
    size_t  get_nb_edges() override;
    size_t  get_nb_vertices() override;
    cpp_int print_num_paths() override;
    void    remove_layers_init();

    PricingSolution farkas_pricing(double* pi) override;
    PricingSolution farkas_pricing(std::span<const double>& pi) override;

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
};

#endif  // PRICER_SOLVER_ZDD_HPP
