//
// Created by daniel on 10/6/21.
//

#include "VariableKeyBase.hpp"
VariableKeyBase::VariableKeyBase(size_t _j, size_t _t, bool _high, bool _root)
    : j(_j),
      t(_t),
      high(_high),
      root(_root) {}

VariableKeyBase::VariableKeyBase() = default;

[[nodiscard]] size_t VariableKeyBase::get_j() const {
    return j;
}
[[nodiscard]] size_t VariableKeyBase::get_t() const {
    return t;
}
[[nodiscard]] bool VariableKeyBase::get_root() const {
    return root;
}
[[nodiscard]] bool VariableKeyBase::get_high() const {
    return high;
}

void VariableKeyBase::set_j(size_t _j) {
    j = _j;
}
void VariableKeyBase::set_t(size_t _t) {
    t = _t;
}
void VariableKeyBase::set_root(bool _root) {
    root = _root;
}
void VariableKeyBase::set_high(bool _high) {
    high = _high;
}

VariableKeyBase::VariableKeyBase(const VariableKeyBase&) = default;
VariableKeyBase::VariableKeyBase(VariableKeyBase&&) noexcept = default;
VariableKeyBase& VariableKeyBase::operator=(const VariableKeyBase&) = default;
VariableKeyBase& VariableKeyBase::operator=(VariableKeyBase&&) noexcept =
    default;

bool VariableKeyBase::operator==(const VariableKeyBase& other) const {
    return (j == other.j && t == other.t && high == other.high);
}
