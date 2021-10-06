//
// Created by daniel on 10/6/21.
//

#include "NodeBdd.hpp"
#include <boost/multiprecision/cpp_int.hpp>  // for cpp_int

void NodeBdd::set_weight(int _weight) {
    weight = _weight;
}
[[nodiscard]] size_t NodeBdd::get_weight() const {
    return weight;
}

[[nodiscard]] Job* NodeBdd::get_job() const {
    return job;
}
[[nodiscard]] size_t NodeBdd::get_nb_job() const {
    return job->job;
}
void NodeBdd::set_job(Job* _job) {
    job = _job;
}

void NodeBdd::add_coeff_list(const std::shared_ptr<BddCoeff>& ptr, bool high) {
    if (high) {
        coeff_list[1].push_back(ptr);
    } else {
        coeff_list[0].push_back(ptr);
    }
}

void NodeBdd::add_coeff_list_clear() {
    for (auto& it : coeff_list) {
        it.clear();
    }
}

std::array<std::vector<NodeBdd::weak_ptr_bddcoeff>, 2>& NodeBdd::get_coeff_list() {
    return coeff_list;
}

bool NodeBdd::operator!=(NodeBdd const& o) const {
    return !(*this == (o));
}

std::ostream& operator<<(std::ostream& os, NodeBdd const& o) {
    os << "(" << o[0];

    for (int i = 1; i < 2; ++i) {
        os << "," << o[i];
    }

    return os << ")";
}

void NodeBdd::init_node(size_t                _weight,
                        [[maybe_unused]] bool _root_node,
                        bool                  _terminal_node) {
    if (!_terminal_node) {
        weight = _weight;
    } else {
        set_job(nullptr);
        weight = -1;
    }
}

void NodeBdd::set_node_id_label(const NodeId& _id) {
    for (auto j = 0UL; j < 2; j++) {
        backward_label[j].set_node_id(_id);
        forward_label[j].set_node_id(_id);
    }
}

void NodeBdd::set_job_label(Job* _job) {
    for (auto j = 0UL; j < 2; j++) {
        backward_label[j].set_job(_job);
        forward_label[j].set_job(_job);
    }
}

auto NodeBdd::operator<=>(const NodeBdd& rhs) {
    return this->forward_label[0].get_f() <=> rhs.forward_label[0].get_f();
}

/** Functions concerning reduced_cost */
NodeBdd::dbl_array& NodeBdd::get_reduced_cost() {
    return reduced_cost;
}
void NodeBdd::reset_reduced_costs() {
    reduced_cost = cost;
}
void NodeBdd::reset_reduced_costs_farkas() {
    reduced_cost = {0.0, 0.0};
}
void NodeBdd::adjust_reduced_costs(double _x, bool high) {
    if (high) {
        reduced_cost[1] -= _x;
    } else {
        reduced_cost[0] -= _x;
    }
}

/** Functions concerning ptr_node_id */
void NodeBdd::set_ptr_node_id(size_t i, size_t j) {
    ptr_node_id = NodeId(i, j);
}
void NodeBdd::set_ptr_node_id(const NodeId& _node_id) {
    ptr_node_id = _node_id;
}
NodeId& NodeBdd::get_ptr_node_id() {
    return ptr_node_id;
}
void NodeBdd::reset_ptr_node_id() {
    ptr_node_id = 0;
}

/** Functions for manipulation of nb_paths */
NodeBdd::big_int& NodeBdd::get_nb_paths() {
    return nb_paths;
}
void NodeBdd::reset_nb_paths() {
    nb_paths = 0UL;
}
void NodeBdd::update_nb_paths(const big_int& x) {
    nb_paths += x;
}

/** Functions for manipulation of visited */
void NodeBdd::reset_visited() {
    visited = false;
}
void NodeBdd::update_visited() {
    visited = true;
}
bool NodeBdd::get_visited() const {
    return visited;
}

/** Functions for manipulation of key */
void NodeBdd::set_key(const size_t& _key) {
    key = _key;
}
size_t& NodeBdd::get_key() {
    return key;
}

/** Functions for manipulation of lp_visited */
void NodeBdd::update_lp_visited(bool _update) {
    lp_visited = _update;
}
bool NodeBdd::get_lp_visited() const {
    return lp_visited;
}

/** Functions for manipulation of lp_x */
[[nodiscard]] NodeBdd::dbl_array & NodeBdd::get_lp_x() {
    return lp_x;
}
double& NodeBdd::get_lp_x(bool _high) {
    return _high ? lp_x[1] : lp_x[0];
}

void NodeBdd::update_lp_x(double _x, bool _high) {
    if (_high) {
        lp_x[1] += _x;
    } else {
        lp_x[0] += _x;
    }
}

void NodeBdd::reset_lp_x() {
    lp_x = {0.0, 0.0};
    lp_visited = false;
}

/** Functions for manipulation of cost */
void NodeBdd::set_cost(double _cost) {
    cost = {0.0, _cost};
}

/** Functions for manipulation of best_sol_x */
void NodeBdd::update_best_sol_x(bool _high) {
    best_sol_x[_high] += 1.0;
}
double NodeBdd::get_best_sol_x(bool _high) {
    return best_sol_x[_high];
}

/** Functions for manipulation of backward distance */
void NodeBdd::reset_backward_distance(int _x) {
    backward_distance = {_x, _x};
}
void NodeBdd::update_backward_distance(
    const std::array<int, 2>& _backward_distance,
    bool                      _high) {
    if (_high) {
        backward_distance[1] =
            ranges::max(_backward_distance) + get_job()->processing_time;
    } else {
        backward_distance[0] = ranges::max(_backward_distance);
    }
}

NodeBdd::int_array& NodeBdd::get_backward_distance() {
    return backward_distance;
}

/** Functions for manipulation of all */
void NodeBdd::reset_all(size_t _nb_elements) {
    all = boost::dynamic_bitset<>{_nb_elements, 0};
}

boost::dynamic_bitset<>& NodeBdd::get_all() {
    return all;
}

void NodeBdd::intersect_all(const boost::dynamic_bitset<>& _set) {
    all &= _set;
}
void NodeBdd::union_all(const boost::dynamic_bitset<>& _set) {
    all |= _set;
}
void NodeBdd::add_element(size_t _element) {
    all[_element] = true;
}
bool NodeBdd::is_element(size_t _element) {
    return all[_element];
}
bool NodeBdd::all_is_empty() {
    return all.empty();
}
size_t NodeBdd::get_first_all() {
    return all.find_first();
}
size_t NodeBdd::get_next_all(size_t _index) {
    return all.find_next(_index);
}

/** Functions for manipulation of calc */
void NodeBdd::reset_calc() {
    calc = {true, true};
}
bool NodeBdd::any_of_calc() {
    return ranges::any_of(calc, [](auto& it) { return it; });
}

bool NodeBdd::get_calc(bool _high) {
    return _high ? calc[1] : calc[0];
}

void NodeBdd::update_calc(bool _high, bool _update) {
    if (_high) {
        calc[1] = _update;
    } else {
        calc[0] = _update;
    }
}

/** Functions for manipulation of in_degree */
int NodeBdd::get_in_degree(bool _high) {
    return _high ? in_degree[1] : in_degree[0];
}

void NodeBdd::update_in_degree(bool _high) {
    if (_high) {
        ++in_degree[1];
    } else {
        ++in_degree[0];
    }
}

bool NodeBdd::alternative_paths() {
    return (ranges::accumulate(in_degree, 0, ranges::plus{}) >= 2);
}

void NodeBdd::reset_in_degree() {
    in_degree = {0, 0};
}

/** Functions for manipulation of Zero-Half cuts formulation cut problem */
void NodeBdd::reset_in_edges() {
    for (auto& it : in_edges) {
        it.clear();
    }
}

std::vector<size_t>& NodeBdd::get_in_edges(bool _high) {
    return _high ? in_edges[1] : in_edges[0];
}
void NodeBdd::add_in_edge(bool _high, size_t _key_parent) {
    if (_high) {
        in_edges[1].emplace_back(_key_parent);
    } else {
        in_edges[0].emplace_back(_key_parent);
    }
}

void NodeBdd::set_key_edge(bool _high, size_t _key) {
    if (_high) {
        key_edges[1] = _key;
    } else {
        key_edges[0] = _key;
    }
}

size_t& NodeBdd::get_key_edge(bool _high) {
    return _high ? key_edges[1] : key_edges[0];
}

void NodeBdd::reset_coeff_cut() {
    coeff_cut = {0.0, 0.0};
}
void NodeBdd::reset_coeff_cut(bool _high) {
    if (_high) {
        coeff_cut[1] = 0.0;
    } else {
        coeff_cut[0] = 0.0;
    }
}

void NodeBdd::update_coeff_cut(bool _high) {
    if (_high) {
        coeff_cut[1] += 1.0;
    } else {
        coeff_cut[0] += 1.0;
    }
}

void NodeBdd::reduce_coeff_cut(bool _high) {
    if (_high) {
        coeff_cut[1] -= 1.0;
    } else {
        coeff_cut[0] -= 1.0;
    }
}

double NodeBdd::get_coeff_cut(bool _high) {
    return _high ? coeff_cut[1] : coeff_cut[0];
}

NodeBdd::grb_array& NodeBdd::get_y() {
    return y;
}
NodeBdd::grb_array& NodeBdd::get_r() {
    return r;
}
void NodeBdd::set_sigma(GRBVar&& _sigma) {
    sigma = _sigma;
}
GRBVar& NodeBdd::get_sigma() {
    return sigma;
}
