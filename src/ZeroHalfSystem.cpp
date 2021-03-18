#include "ZeroHalfSystem.hpp"
#include <algorithm>
#include <boost/container_hash/hash_fwd.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstddef>
#include <limits>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/min.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/stride.hpp>
#include <range/v3/view/transform.hpp>
#include <unordered_map>
#include <vector>

namespace vs = ranges::views;

ZeroHalfSystem::ZeroHalfSystem(const MatrixDouble& _A,
                               const VectorDouble& _b,
                               const VectorDouble& _x)
    : x_star(_x),
      row_index(_A.size(), boost::dynamic_bitset{_A.size()}) {
    auto min_elem =
        ranges::min(ranges::min(ranges::views::transform(
                        _A, [](auto& row) { return ranges::min(row); })),
                    ranges::min(_b));

    auto add_elem = 0.0;
    if (min_elem < 0.0) {
        add_elem = std::ceil(-min_elem / HALF) * HALF;
    }

    auto func_add_value = [&add_elem](auto& it) -> int {
        return static_cast<int>(add_elem + it) % 2;
    };

    A_bar = vs::transform(_A,
                          [&func_add_value](auto& row) {
                              return vs::transform(row, func_add_value) |
                                     ranges::to<std::vector>();
                          }) |
            ranges::to<std::vector>();

    b_bar = vs::transform(_b, func_add_value) | ranges::to<std::vector>;

    auto tmp = vs::transform(
        _A, [&_x](auto& row) { return ranges::inner_product(row, _x, 0.0); });

    slack =
        vs::zip_with([](auto lhs, auto rhs) { return lhs - rhs; }, _b, tmp) |
        ranges::to<std::vector>();

    for (auto&& [i, set] : row_index | vs::enumerate) {
        set[i] = true;
    }

    reduce_system();
}

void ZeroHalfSystem::remove_row(size_t _row) {
    ranges::swap(A_bar[_row], A_bar[nb_rows - 1]);
    ranges::swap(slack[_row], slack[nb_rows - 1]);
    ranges::swap(b_bar[_row], b_bar[nb_rows - 1]);
    ranges::swap(row_index[_row], row_index[nb_rows - 1]);
    --nb_rows;
    // A_bar.pop_back();
}

void ZeroHalfSystem::remove_col(size_t _col) {
    for (auto& it : A_bar) {
        ranges::swap(it[_col], it[nb_columns - 1]);
    }
    ranges::swap(x_star[_col], x_star[nb_columns - 1]);
    --nb_columns;
}

void ZeroHalfSystem::reduce_system() {
    auto col = 0UL;

    /** Rule 1 */
    while (col < nb_columns) {
        if (x_star[col] < EPS) {
            remove_col(col);
        } else {
            ++col;
        }
    }

    /** Rule 2 */
    auto row = 0UL;
    while (row < nb_rows) {
        auto any_one =
            ranges::any_of(A_bar[row], [](auto& a) { return a == 1; }) ||
            (b_bar[row] == 1);

        if (any_one) {
            ++row;
        } else {
            remove_row(row);
        }
    }

    /** Rule 3,4,5 */

    std::unordered_map<size_t, int> column_dict{};
    col = 0UL;
    while (col < nb_columns) {
        size_t key = 0UL;
        int    counter{};
        int    last_one{-1};
        for (auto&& [r, a] : A_bar | vs::enumerate) {
            if (a[col]) {
                boost::hash_combine(key, r);
                ++counter;
                last_one = r;
            }
        }

        if (counter == 0) {
            remove_col(col);
        } else if (counter == 1) {
            slack[last_one] += x_star[col];
            remove_col(col);
        } else if (column_dict.find(key) != column_dict.end()) {
            auto aux = column_dict[key];
            x_star[aux] += x_star[col];
            remove_col(col);
        } else {
            column_dict[key] = col;
            ++col;
        }
    }

    /** Rule 6 */
    row = 0UL;
    while (row < nb_rows) {
        if (slack[row] >= 1.0) {
            remove_row(row);
        } else {
            ++row;
        }
    }
    /** Use hash combine to detect equivalent columns */
    /** Rule 7 Dekoster et al. */
    std::unordered_map<size_t, int> row_dict;
    row = 0UL;
    while (row < nb_rows) {
        size_t key = 0UL;
        ranges::for_each(A_bar[row], [&key](auto& it_col) {
            if (it_col) {
                boost::hash_combine(key, it_col);
            }
        });

        if (b_bar[row]) {
            boost::hash_combine(key, nb_columns);
        }

        if (row_dict.find(key) == row_dict.end()) {
            row_dict[key] = row;
            ++row;
        } else {
            auto tmp = row_dict[key];
            if (slack[tmp] < slack[row]) {
                remove_row(row);
            } else {
                ranges::swap(row_index[tmp], row_index[row]);
                ranges::swap(slack[tmp], slack[row]);
                remove_row(row);
            }
        }
    }
}

void ZeroHalfSystem::reduce_gauss() {
    auto row = 0UL;
    while (row < nb_rows) {
        if (slack[row] < EPS) {
            auto high_x = std::numeric_limits<double>::min();
            auto best_col = -1;
            auto counter = 0;
            for (auto j = 0UL; j < nb_columns; ++j) {
                if (A_bar[row][j]) {
                    ++counter;
                    if (x_star[j] > high_x) {
                        high_x = x_star[j];
                        best_col = j;
                    }
                }
            }

            if (counter == 0) {
                if (b_bar[row]) {
                    /** construct ineq */
                }
                remove_row(row);
            } else {
                for (auto&& [r, x] : vs::enumerate(A_bar)) {
                    if (r != row && x[best_col] == 1) {
                        add_to_row(row, r);
                    }
                    slack[row] += x_star[best_col];
                    remove_col(best_col);
                    ++row;
                }
            }
        } else {
            ++row;
        }
    }
}

void ZeroHalfSystem::evaluate_rows(const std::vector<int>& _rows) {
    auto sum_slack = 0.0;
    ranges::for_each(_rows, [&](int idx) { sum_slack += slack[idx]; });
    if (sum_slack <= 1.0 - HALF * EPS) {
        std::vector<int> v(nb_rows, 0);
        ranges::for_each(_rows, [&](int idx) { v[idx] = 1; });
        auto odd = (ranges::inner_product(v, b_bar, 0) % 2) == 1;

        // auto test = A_bar | vs::join | ;
        if (odd) {
            auto left = vs::ints(0UL, nb_columns) | vs::transform([&](int i) {
                            return ranges::inner_product(
                                       v | vs::all,
                                       A_bar | vs::join | vs::drop(i) |
                                           vs::stride(nb_columns),
                                       0) %
                                   2;
                        });
            auto val = ranges::inner_product(left, x_star, 0.0);
            if (val < 1 - HALF * EPS) {
                auto result = row_index[0];
                for (auto& it : row_index | vs::drop(0)) {
                    result ^= it;
                }
                /** construct ineq */
            }
        }
    }
}

void ZeroHalfSystem::add_to_row(size_t i, size_t j) {
    for (auto&& [lhs, rhs] : vs::zip(A_bar[i], A_bar[j])) {
        rhs = (lhs + rhs) % 2;
    }
    b_bar[j] = (b_bar[i] + b_bar[j]) % 2;
    row_index[j] = row_index[i] ^ row_index[j];
    slack[j] = slack[i] + slack[j];
}
