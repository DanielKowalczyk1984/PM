//
// Created by daniel on 10/6/21.
//
#include <cstddef>
#include <ostream>
#include <boost/container_hash/extensions.hpp>  // for hash_combine
#include "VariableKeyBase.hpp"
#ifndef PM_BDDCOEFF_H
#define PM_BDDCOEFF_H

class BddCoeff : public VariableKeyBase {
   private:
    size_t row;
    double coeff;
    double value;

   public:
    BddCoeff(size_t _j,
             size_t _t,
             double _coeff,
             double _value = 0.0,
             size_t _row = 0UL,
             bool   _high = true,
             bool   _root = false);

    BddCoeff(const BddCoeff&);
    BddCoeff& operator=(const BddCoeff&);
    BddCoeff(BddCoeff&& op) noexcept ;
    BddCoeff& operator=(BddCoeff&& op) noexcept ;
    ~BddCoeff() override = default;

    [[nodiscard]] double get_coeff() const;
    [[nodiscard]] double get_value() const;
    [[nodiscard]] size_t get_row() const;
    void set_value(double _value);
    void set_row(size_t _row);

    bool operator==(const BddCoeff& other);

    friend std::ostream& operator<<(std::ostream& os, const BddCoeff& object);
    friend bool operator==(const BddCoeff& lhs, const BddCoeff& rhs);
};

namespace std {
template <>
struct hash<BddCoeff> {
    std::size_t operator()(auto const& s) const noexcept {
        std::size_t seed = 0;
        boost::hash_combine(seed, s.get_j());
        boost::hash_combine(seed, s.get_t());
        boost::hash_combine(seed, s.get_high());

        return seed;  // or use boost::hash_combine
    }
};
}  // namespace std
#endif  // PM_BDDCOEFF_H
