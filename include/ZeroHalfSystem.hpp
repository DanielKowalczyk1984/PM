#ifndef __ZEROHALFSYSTEM_H__
#define __ZEROHALFSYSTEM_H__

#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstddef>
#include <vector>
class ZeroHalfSystem {
    using MatrixDouble = std::vector<std::vector<double>>;
    using VectorDouble = std::vector<double>;

   public:
    ZeroHalfSystem() = default;
    ZeroHalfSystem(const MatrixDouble& _A,
                   const VectorDouble& b,
                   const VectorDouble& _x);
    ZeroHalfSystem(ZeroHalfSystem&&) = default;
    ZeroHalfSystem(const ZeroHalfSystem&) = default;
    ZeroHalfSystem& operator=(ZeroHalfSystem&&) = default;
    ZeroHalfSystem& operator=(const ZeroHalfSystem&) = default;
    ~ZeroHalfSystem() = default;

   private:
    std::vector<std::vector<int>> A_bar;
    std::vector<int>              b_bar;
    std::vector<double>           x_star;
    std::vector<double>           slack;

    std::vector<boost::dynamic_bitset<>> row_index;
    size_t                               nb_rows{};
    size_t                               nb_columns{};

    static constexpr double HALF = 2.0;
    static constexpr double EPS = 1e-6;

    void remove_row(size_t _row);
    void remove_col(size_t _col);

    void reduce_system();
    void reduce_gauss();
    void evaluate_rows(const std::vector<int> _rows);

    void add_to_row(size_t i, size_t j);
};

#endif  // __ZEROHALFSYSTEM_H__