//
// Created by daniel on 10/6/21.
//

#include "BddCoeff.hpp"

BddCoeff::BddCoeff(size_t _j,
                   size_t _t,
                   double _coeff,
                   double _value,
                   size_t _row,
                   bool   _high,
                   bool   _root)
    : VariableKeyBase(_j, _t, _high, _root),
      row(_row),
      coeff(_coeff),
      value(_value) {}

BddCoeff::BddCoeff(const BddCoeff&) = default;
BddCoeff& BddCoeff::operator=(const BddCoeff&) = default;
BddCoeff::BddCoeff(BddCoeff&& op) noexcept = default;
BddCoeff& BddCoeff::operator=(BddCoeff&& op) noexcept = default;

std::ostream& operator<<(std::ostream& os, const BddCoeff& object) {
    return os << "(j = " << object.get_j() << ", t = " << object.get_t()
              << ", x = " << object.get_value()
              << ", high = " << object.get_high() << " )\n";
}

bool operator==(const BddCoeff& lhs, const BddCoeff& rhs) {
    return lhs.get_j() == rhs.get_j() && lhs.get_t() == rhs.get_t() &&
           lhs.get_high() == rhs.get_high();
}

bool BddCoeff::operator==(const BddCoeff& other) {
    return get_j() == other.get_j() && get_t() == other.get_t() &&
           get_high() == other.get_high();
}

[[nodiscard]] double BddCoeff::get_coeff() const {
    return coeff;
}

[[nodiscard]] double BddCoeff::get_value() const {
    return value;
}

void BddCoeff::set_value(double _value) {
    value = _value;
}

void BddCoeff::set_row(size_t _row) {
    row = _row;
}

[[nodiscard]] size_t BddCoeff::get_row() const {
    return row;
}
