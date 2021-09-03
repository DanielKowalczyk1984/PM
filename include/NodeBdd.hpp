#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP

#include <gurobi_c++.h>                             // for GRBVar
#include <array>                                    // for array
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/multiprecision/cpp_int.hpp>         // for cpp_int
#include <cstddef>                                  // for size_t
#include <memory>                                   // for allocator, shared...
#include <ostream>                                  // for operator<<, ostream
#include <range/v3/algorithm/any_of.hpp>            // for any_of
#include <range/v3/algorithm/max.hpp>               // for max, max_fn
#include <range/v3/functional/arithmetic.hpp>       // for plus
#include <range/v3/functional/comparisons.hpp>      // for less
#include <range/v3/functional/identity.hpp>         // for identity
#include <range/v3/numeric/accumulate.hpp>          // for accumulate, accum...
#include <vector>                                   // for vector
#include "Job.h"                                    // for Job
#include "Label.hpp"                                // for Label
#include "NodeBase.hpp"                             // for NodeBase
#include "NodeId.hpp"                               // for NodeId

class BddCoeff;

template <typename T = double>
class NodeBdd : public NodeBase {
    using dbl_array = std::array<double, 2>;
    using int_array = std::array<int, 2>;
    using bool_array = std::array<bool, 2>;
    using big_int = boost::multiprecision::cpp_int;
    using weak_ptr_nodeid = std::weak_ptr<NodeId>;
    using weak_ptr_bddcoeff = std::weak_ptr<BddCoeff>;

   private:
    /** Job chatacteristics */
    Job*   job{nullptr};
    int    weight{};
    size_t key{};

    bool   visited{false};
    NodeId ptr_node_id{};

    /** Calculation jobs */
    boost::dynamic_bitset<> all{};

    /** original model data */
    dbl_array lp_x{0.0, 0.0};
    dbl_array cost{0.0, 0.0};
    dbl_array best_sol_x{0.0, 0.0};
    dbl_array reduced_cost{0.0, 0.0};

    std::array<std::vector<weak_ptr_bddcoeff>, 2> coeff_list{};

    /** Evaluation of backward distance */
    int_array backward_distance{};

    /** Evaluation nb_paths */
    big_int nb_paths{};

    /** Reducing size of decision diagram */
    bool_array calc{true, true};

    /** Searching equivalent paths */
    int_array in_degree{};

    /** Zero-Half cuts formulation variables */
    bool                                        lp_visited{false};
    std::array<std::vector<weak_ptr_nodeid>, 2> in_edges{};
    dbl_array                                   coeff_cut{0.0, 0.0};
    std::array<GRBVar, 2>                       y{};
    std::array<GRBVar, 2>                       r{};
    GRBVar                                      sigma{};

   public:
    std::array<Label<NodeBdd<T>, T>, 2> forward_label{};
    std::array<Label<NodeBdd<T>, T>, 2> backward_label{};

    /**
     * Constructor
     */
    NodeBdd() = default;
    NodeBdd(size_t i, size_t j) : NodeBase(i, j) {}

    NodeBdd(const NodeBdd<T>& src) = default;
    NodeBdd(NodeBdd<T>&& src) noexcept = default;
    NodeBdd<T>& operator=(const NodeBdd<T>& src) = default;
    NodeBdd<T>& operator=(NodeBdd<T>&& src) noexcept = default;
    ~NodeBdd() = default;

    void              set_weight(int _weight) { weight = _weight; }
    [[nodiscard]] int get_weight() const { return weight; }

    [[nodiscard]] Job*          get_job() const { return job; }
    [[nodiscard]] inline size_t get_nb_job() const { return job->job; }
    void                        set_job(Job* _job) { job = _job; }

    void add_coeff_list(const std::shared_ptr<BddCoeff>& ptr, bool high) {
        if (high) {
            coeff_list[1].push_back(ptr);
        } else {
            coeff_list[0].push_back(ptr);
        }
    }

    void add_coeff_list_clear() {
        for (auto& it : coeff_list) {
            it.clear();
        }
    }

    auto& get_coeff_list() { return coeff_list; }

    bool operator!=(NodeBdd const& o) const { return !operator==(o); }

    friend std::ostream& operator<<(std::ostream& os, NodeBdd const& o) {
        os << "(" << o[0];

        for (int i = 1; i < 2; ++i) {
            os << "," << o[i];
        }

        return os << ")";
    }

    void init_node(int                   _weight,
                   [[maybe_unused]] bool _root_node = false,
                   bool                  _terminal_node = false) {
        if (!_terminal_node) {
            weight = _weight;
        } else {
            set_job(nullptr);
            weight = -1;
        }
    }

    void set_node_id_label(const NodeId& _id) {
        for (auto j = 0UL; j < 2; j++) {
            backward_label[j].set_node_id(_id);
            forward_label[j].set_node_id(_id);
        }
    }

    void set_job_label(Job* _job) {
        for (auto j = 0UL; j < 2; j++) {
            backward_label[j].set_job(_job);
            forward_label[j].set_job(_job);
        }
    }

    friend bool operator<=>(const NodeBdd<T>& lhs, const NodeBdd<T>& rhs) {
        return lhs.forward_label[0] <=> rhs.forward_label[0].f;
    }

    /** Functions concerning reduced_cost */
    inline auto& get_reduced_cost() { return reduced_cost; }
    void         reset_reduced_costs() { reduced_cost = cost; }
    void         reset_reduced_costs_farkas() { reduced_cost = {0.0, 0.0}; }
    void         adjust_reduced_costs(double _x, bool high) {
        if (high) {
            reduced_cost[1] -= _x;
        } else {
            reduced_cost[0] -= _x;
        }
    }

    /** Functions concerning ptr_node_id */
    void  set_ptr_node_id(size_t i, size_t j) { ptr_node_id = NodeId(i, j); }
    void  set_ptr_node_id(const NodeId& _node_id) { ptr_node_id = _node_id; }
    auto& get_ptr_node_id() { return ptr_node_id; }
    void  reset_ptr_node_id() { ptr_node_id = 0; }

    /** Functions for manipulation of nb_paths */
    big_int& get_nb_paths() { return nb_paths; }
    void     reset_nb_paths() { nb_paths = 0UL; }
    void     update_nb_paths(const big_int& x = 1) { nb_paths += x; }

    /** Functions for manipulation of visited */
    void reset_visited() { visited = false; }
    void update_visited() { visited = true; }
    bool get_visited() { return visited; }

    /** Functions for manipulation of key */
    void    set_key(const size_t& _key) { key = _key; }
    size_t& get_key() { return key; }

    /** Functions for manipulation of lp_visited */
    void update_lp_visited(bool _update) { lp_visited = _update; }
    bool get_lp_visited() { return lp_visited; }

    /** Functions for manipulation of lp_x */
    [[nodiscard]] dbl_array& get_lp_x() { return lp_x; }

    void update_lp_x(double _x, bool _high) {
        if (_high) {
            lp_x[1] += _x;
        } else {
            lp_x[0] += _x;
        }
    }

    void reset_lp_x() {
        lp_x = {0.0, 0.0};
        lp_visited = false;
    }

    /** Functions for manipulation of cost */
    void set_cost(double _cost) { cost = {0.0, _cost}; }

    /** Functions for manipulation of best_sol_x */
    void   update_best_sol_x(bool _high) { best_sol_x[_high] += 1.0; }
    double get_best_sol_x(bool _high) { return best_sol_x[_high]; }

    /** Functions for manipulation of backward distance */
    void reset_backward_distance(int _x = 0) { backward_distance = {_x, _x}; }
    void update_backward_distance(const std::array<int, 2>& _backward_distance,
                                  bool                      _high) {
        if (_high) {
            backward_distance[1] =
                ranges::max(_backward_distance) + get_job()->processing_time;
        } else {
            backward_distance[0] = ranges::max(_backward_distance);
        }
    }

    std::array<int, 2>& get_backward_distance() { return backward_distance; }

    /** Functions for manipulation of all */
    void reset_all(size_t _nb_elements) {
        all = boost::dynamic_bitset<>{_nb_elements, 0};
    }

    boost::dynamic_bitset<>& get_all() { return all; }

    void   intersect_all(const boost::dynamic_bitset<>& _set) { all &= _set; }
    void   union_all(const boost::dynamic_bitset<>& _set) { all |= _set; }
    void   add_element(size_t _element) { all[_element] = true; }
    bool   is_element(size_t _element) { return all[_element]; }
    bool   all_is_empty() { return all.empty(); }
    size_t get_first_all() { return all.find_first(); }
    size_t get_next_all(size_t _index) { return all.find_next(_index); }

    /** Functions for manipulation of calc */
    void reset_calc() { calc = {true, true}; }
    bool any_of_calc() {
        return ranges::any_of(calc, [](auto& it) { return it; });
    }

    bool get_calc(bool _high) { return _high ? calc[1] : calc[0]; }

    void update_calc(bool _high, bool _update = false) {
        if (_high) {
            calc[1] = _update;
        } else {
            calc[0] = _update;
        }
    }

    /** Functions for manipulation of in_degree */
    int get_in_degree(bool _high) {
        return _high ? in_degree[1] : in_degree[0];
    }

    void update_in_degree(bool _high) {
        if (_high) {
            ++in_degree[1];
        } else {
            ++in_degree[0];
        }
    }

    bool alternative_paths() {
        return (ranges::accumulate(in_degree, 0, ranges::plus{}) >= 2);
    }

    void reset_in_degree() { in_degree = {0, 0}; }

    /** Functions for manipulation of Zero-Half cuts formulation cut problem */
    void reset_in_edges() {
        for (auto& it : in_edges) {
            it.clear();
        }
    }

    auto& get_in_edges(bool _high) { return _high ? in_edges[1] : in_edges[0]; }

    void reset_coeff_cut() { coeff_cut = {0.0, 0.0}; }
    void reset_coeff_cut(bool _high) {
        if (_high) {
            coeff_cut[1] = 0.0;
        } else {
            coeff_cut[0] = 0.0;
        }
    }

    void update_coeff_cut(bool _high) {
        if (_high) {
            coeff_cut[1] += 1.0;
        } else {
            coeff_cut[0] += 1.0;
        }
    }

    void reduce_coeff_cut(bool _high) {
        if (_high) {
            coeff_cut[1] -= 1.0;
        } else {
            coeff_cut[0] -= 1.0;
        }
    }

    double get_coeff_cut(bool _high) {
        return _high ? coeff_cut[1] : coeff_cut[0];
    }

    auto& get_y() { return y; }
    auto& get_r() { return r; }
    void  set_sigma(GRBVar&& _sigma) { sigma = _sigma; }
    auto& get_sigma() { return sigma; }
};

#endif  // NODE_DURATION_HPP
