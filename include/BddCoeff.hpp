//
// Created by daniel on 10/6/21.
//
#include "VariableKeyBase.h"
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
             bool   _root = false)
        : VariableKeyBase(_j, _t, _high, _root),
          row(_row),
          coeff(_coeff),
          value(_value){};

    BddCoeff(const BddCoeff&) = default;
    BddCoeff& operator=(const BddCoeff&) = default;
    BddCoeff(BddCoeff&& op) = default;
    BddCoeff& operator=(BddCoeff&& op) = default;
    ~BddCoeff() override = default;

    [[nodiscard]] inline double get_coeff() const { return coeff; }

    [[nodiscard]] inline double get_value() const { return value; }

    inline void set_value(double _value) { value = _value; }

    inline void set_row(size_t _row) { row = _row; }

    [[nodiscard]] inline size_t get_row() const { return row; }

    friend bool operator==(const BddCoeff& lhs, const BddCoeff& rhs) {
        return lhs.get_j() == rhs.get_j() && lhs.get_t() == rhs.get_t() &&
               lhs.get_high() == rhs.get_high();
    };

    bool operator==(const BddCoeff& other) {
        return get_j() == other.get_j() && get_t() == other.get_t() &&
               get_high() == other.get_high();
    }

    friend std::ostream& operator<<(std::ostream& os, const BddCoeff& object) {
        return os << "(j = " << object.get_j() << ", t = " << object.get_t()
                  << ", x = " << object.get_value()
                  << ", high = " << object.get_high() << " )\n";
    }
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
}
#endif  // PM_BDDCOEFF_H
