#ifndef _MODEL_INTERFACE
#define _MODEL_INTERFACE

#include <bits/c++config.h>
#include <fmt/core.h>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <unordered_map>
#include <vector>
namespace vs = ranges::views;

class VariableKeyBase {
   private:
    int  j{-1};
    int  t{-1};
    bool high{false};
    bool root{false};

   public:
    VariableKeyBase(int _j, int _t, bool _high = true, bool _root = false)
        : j(_j),
          t(_t),
          high(_high),
          root(_root) {}

    VariableKeyBase() = default;

    [[nodiscard]] inline int get_j() const { return j; }

    [[nodiscard]] inline int get_t() const { return t; }

    [[nodiscard]] inline bool get_root() const { return root; }

    [[nodiscard]] inline bool get_high() const { return high; }

    inline void set_j(int _j) { j = _j; }

    inline void set_t(int _t) { t = _t; }

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

class ConstraintBase {
    char   sense;
    double rhs;
    bool   can_be_deleted;

   public:
    [[nodiscard]] inline double get_rhs() const { return rhs; }

    [[nodiscard]] inline char get_sense() const { return sense; }

    [[nodiscard]] inline bool get_can_be_deleted() const {
        return can_be_deleted;
    }

    ConstraintBase(const ConstraintBase&) = default;
    ConstraintBase(ConstraintBase&&) = default;
    ConstraintBase& operator=(const ConstraintBase&) = default;
    ConstraintBase& operator=(ConstraintBase&&) = default;
    virtual ~ConstraintBase() = default;

    ConstraintBase(char _sense, double _rhs, bool _can_be_delete = false)
        : sense(_sense),
          rhs(_rhs),
          can_be_deleted(_can_be_delete) {}

    virtual double operator()(const VariableKeyBase&) = 0;
};

class ConstraintAssignment : public ConstraintBase {
   private:
    int row;

   public:
    explicit ConstraintAssignment(int _row)
        : ConstraintBase('>', 1.0),
          row(_row) {}

    double operator()(const VariableKeyBase& key) override {
        if (key.get_j() == row && key.get_high()) {
            return 1.0;
        }

        return 0.0;
    }

    ConstraintAssignment(const ConstraintAssignment&) = default;
    ConstraintAssignment(ConstraintAssignment&&) = default;
    ConstraintAssignment& operator=(const ConstraintAssignment&) = default;
    ConstraintAssignment& operator=(ConstraintAssignment&&) = default;
    ~ConstraintAssignment() override = default;
};

class ConstraintConvex : public ConstraintBase {
   public:
    explicit ConstraintConvex(double _rhs) : ConstraintBase('>', _rhs) {}

    double operator()(const VariableKeyBase& key) override {
        if (!key.get_t()) {
            return -1.0;
        }
        return 0.0;
    }
};

class ReformulationModel : public std::vector<std::shared_ptr<ConstraintBase>> {
   public:
    ReformulationModel(int nb_assignments, int nb_machines);
    ReformulationModel(ReformulationModel&&) = default;
    ReformulationModel& operator=(ReformulationModel&&) = default;
    ReformulationModel(const ReformulationModel&) = default;
    ReformulationModel& operator=(const ReformulationModel&) = default;
    ~ReformulationModel() = default;

    inline void delete_constraint(auto c) {
        if ((*this)[c]->get_can_be_deleted()) {
            (*this)[c].reset();
        }
    }

    inline void delete_constraints(int first, int nb_del) {
        auto it = this->begin() + first;
        this->erase(it, it + nb_del);
    }
};

class BddCoeff : public VariableKeyBase {
   private:
    size_t row;
    double coeff;
    double value;

   public:
    BddCoeff(int    _j,
             int    _t,
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

class GenericData : public std::unordered_map<VariableKeyBase, double> {
   private:
    static constexpr double EPS_GENERIC_DATA = 1e-6;

   public:
    GenericData() = default;
    ~GenericData() = default;
    GenericData(GenericData&&) = default;  // movable and noncopyable
    GenericData(const GenericData&) = default;
    GenericData& operator=(GenericData&&) = default;
    GenericData& operator=(const GenericData&) = default;

    void add_coeff_hash_table(int _j, int _t, bool _high, double _coeff) {
        VariableKeyBase key(_j, _t, _high);

        auto it = this->find(key);
        if (it == this->end()) {
            (*this)[key] = _coeff;
        } else {
            (*this)[key] += _coeff;
        }
    }

    void list_coeff() {
        for (auto& it : (*this)) {
            fmt::print("{} ({},{},{})", it.second, it.first.get_j(),
                       it.first.get_t(), it.first.get_high());
        }
        fmt::print("\n");
    }

    friend bool operator==(const GenericData& lhs, const GenericData& rhs) {
        if (lhs.size() != rhs.size()) {
            return false;
        }

        for (const auto& it1 : lhs) {
            const auto it2 = rhs.find(it1.first);
            if (it2 == rhs.end()) {
                return false;
            }

            if (fabs(it1.second - (*it2).second) > EPS_GENERIC_DATA) {
                return false;
            }
        }

        return true;
    };

    friend bool operator!=(const GenericData& lhs, const GenericData& rhs) {
        return !(lhs == rhs);
    }
};

class ConstraintGeneric : public ConstraintBase {
   private:
    std::shared_ptr<GenericData> data;

   public:
    ConstraintGeneric(GenericData* _data,
                      double       _rhs,
                      char         _sense = '>',
                      bool         _can_be_deleted = true)
        : ConstraintBase(_sense, _rhs, _can_be_deleted),
          data(_data) {}

    explicit ConstraintGeneric(double _rhs,
                               char   _sense = '>',
                               bool   _can_be_deleted = true)
        : ConstraintBase(_sense, _rhs, _can_be_deleted),
          data(nullptr) {}

    ~ConstraintGeneric() override = default;
    ConstraintGeneric(ConstraintGeneric&& op) = default;
    ConstraintGeneric& operator=(ConstraintGeneric&& op) = default;
    ConstraintGeneric& operator=(const ConstraintGeneric&) = default;
    ConstraintGeneric(const ConstraintGeneric&) = default;

    double operator()(const VariableKeyBase& key) override {
        auto it = data->find(key);
        if (it == data->end()) {
            return 0.0;
        } else {
            return (*it).second;
        }
    }

    friend bool operator!=(const ConstraintGeneric& lhs,
                           const ConstraintGeneric& rhs) {
        return !(*lhs.data == *rhs.data);
    }

    friend bool operator==(const ConstraintGeneric& lhs,
                           const ConstraintGeneric& rhs) {
        return (*lhs.data == *rhs.data);
    }

    void list_coeff() { data->list_coeff(); }
};

template <typename T = BddCoeff>
class OriginalConstraint {
   private:
    size_t                        id_constr{};
    std::weak_ptr<ConstraintBase> constr;
    std::list<std::shared_ptr<T>> coeff_list;

   public:
    explicit OriginalConstraint(const std::shared_ptr<ConstraintBase>& _constr)
        : constr(_constr){};
    OriginalConstraint() : constr(){};
    ~OriginalConstraint() = default;
    OriginalConstraint(OriginalConstraint&&) noexcept = default;
    OriginalConstraint& operator=(OriginalConstraint&&) noexcept = default;
    OriginalConstraint<T>(const OriginalConstraint<T>&) = default;
    OriginalConstraint<T>& operator=(const OriginalConstraint<T>&) = default;

    inline std::list<std::shared_ptr<T>>& get_coeff_list() {
        return coeff_list;
    }

    inline ConstraintBase* get_constr() {
        auto aux = constr.lock();
        if (aux) {
            return aux.get();
        } else {
            return nullptr;
        }
    }

    inline void add_coeff_to_list(std::shared_ptr<T> _coeff) {
        coeff_list.push_back(_coeff);
    }

    void clear_coeff() { coeff_list.clear(); }
};

template <typename T = BddCoeff>
class OriginalModel : public std::vector<OriginalConstraint<T>> {
   public:
    explicit OriginalModel(const ReformulationModel& model)
        : std::vector<OriginalConstraint<T>>(
              vs::iota(0UL, model.size()) | vs::transform([&](auto i) {
                  return OriginalConstraint<T>(model[i]);
              }) |
              ranges::to<std::vector<OriginalConstraint<T>>>()) {}

    OriginalModel(OriginalModel&& op) noexcept = default;
    OriginalModel& operator=(OriginalModel&& op) noexcept = default;
    OriginalModel(const OriginalModel&) = default;
    OriginalModel& operator=(const OriginalModel&) = default;
    ~OriginalModel() = default;

    void add_coeff_list(int c, std::shared_ptr<T> coeff) {
        (*this)[c].add_coeff_to_list(coeff);
    }

    ConstraintBase* get_constraint(int c) { return (*this)[c].get_constr(); }

    inline std::list<std::shared_ptr<T>>& get_coeff_list(int c) {
        return (*this)[c].get_coeff_list();
    }

    inline void add_constraint(const std::shared_ptr<ConstraintBase>& _constr) {
        this->push_back(OriginalConstraint<>(_constr));
    }

    inline size_t get_nb_constraints() { return this->size(); }

    inline void delete_constraints(int first, int nb_del) {
        auto it = this->begin() + first;
        this->erase(it, it + nb_del);
    }

    void clear_all_coeff() {
        for (auto& it : *this) {
            it.clear_coeff();
        }
    }
};

#endif