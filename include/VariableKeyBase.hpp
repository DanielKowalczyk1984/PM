//
// Created by daniel on 10/6/21.
//

#ifndef PM_VARIABLEKEYBASE_H
#define PM_VARIABLEKEYBASE_H

class VariableKeyBase {
   private:
    size_t j{};
    size_t t{};
    bool   high{false};
    bool   root{false};

   public:
    VariableKeyBase(size_t _j, size_t _t, bool _high = true, bool _root = false)
        : j(_j),
          t(_t),
          high(_high),
          root(_root) {}

    VariableKeyBase() = default;

    [[nodiscard]] inline size_t get_j() const { return j; }

    [[nodiscard]] inline size_t get_t() const { return t; }

    [[nodiscard]] inline bool get_root() const { return root; }

    [[nodiscard]] inline bool get_high() const { return high; }

    inline void set_j(size_t _j) { j = _j; }

    inline void set_t(size_t _t) { t = _t; }

    inline void set_root(bool _root) { root = _root; }

    inline void set_high(bool _high) { high = _high; }

    VariableKeyBase(const VariableKeyBase&) = default;
    VariableKeyBase(VariableKeyBase&&) = default;
    VariableKeyBase& operator=(const VariableKeyBase&) = default;
    VariableKeyBase& operator=(VariableKeyBase&&) = default;
    virtual ~VariableKeyBase() = default;

    bool operator==(const VariableKeyBase& other) const {
        return (j == other.j && t == other.t && high == other.high);
    }

    // friend bool operator==(const VariableKeyBase& lhs,
    //                        const VariableKeyBase& rhs) {
    //     return (lhs.j == rhs.j && lhs.t == rhs.t && lhs.high == rhs.high);
    // }
};

namespace std {
template <>
struct hash<VariableKeyBase> {
    std::size_t operator()(auto const& s) const noexcept {
        std::size_t seed = 0;
        boost::hash_combine(seed, s.get_j());
        boost::hash_combine(seed, s.get_t());
        boost::hash_combine(seed, s.get_high());

        return seed;  // or use boost::hash_combine
    }
};
}  // namespace std
#endif  // PM_VARIABLEKEYBASE_H
