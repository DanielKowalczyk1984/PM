#ifndef _MODEL_INTERFACE
#define _MODEL_INTERFACE

#include <bits/c++config.h>
#include <fmt/core.h>
#include <NodeId.hpp>
#include <boost/container_hash/hash_fwd.hpp>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include "wctparms.h"
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

    inline bool get_root() { return root; }

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
};

class ConstraintBase {
    double rhs;
    char   sense;
    bool   can_be_deleted;

   public:
    inline double get_rhs() { return rhs; }

    inline char get_sense() { return sense; }

    inline bool get_can_be_deleted() { return can_be_deleted; }

    ConstraintBase(const ConstraintBase&) = default;
    ConstraintBase(ConstraintBase&&) = default;
    ConstraintBase& operator=(const ConstraintBase&) = default;
    ConstraintBase& operator=(ConstraintBase&&) = default;
    virtual ~ConstraintBase() = default;

    ConstraintBase(char _sense, double _rhs, bool _can_be_delete = false)
        : sense(_sense),
          rhs(_rhs),
          can_be_deleted(_can_be_delete) {}

    virtual double get_var_coeff(VariableKeyBase*) = 0;
};

class ConstraintAssignment : public ConstraintBase {
   private:
    int job;

   public:
    explicit ConstraintAssignment(int _job)
        : ConstraintBase('>', 1.0),
          job(_job) {}

    double get_var_coeff(VariableKeyBase* key) override {
        if (key->get_j() == job && key->get_high()) {
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

    double get_var_coeff(VariableKeyBase* key) override {
        if (key->get_t() == 0) {
            return -1.0;
        }

        return 0.0;
    }
};

class ReformulationModel {
   private:
    std::vector<std::shared_ptr<ConstraintBase>> constraint_array;

   public:
    ReformulationModel(int nb_assignments, int nb_machines);
    ~ReformulationModel() = default;
    ReformulationModel(ReformulationModel&&) noexcept =
        default;  // movable and noncopyable ReformulationModel&
    ReformulationModel& operator=(ReformulationModel&&) = default;
    ReformulationModel(const ReformulationModel&) = default;
    ReformulationModel& operator=(const ReformulationModel&) = default;

    [[nodiscard]] inline size_t get_nb_constraints() const {
        return constraint_array.size();
    };

    [[nodiscard]] inline ConstraintBase* get_constraint(int c) const {
        return constraint_array[c].get();
    };

    inline std::shared_ptr<ConstraintBase> get_constraint_ptr(int c) const {
        return constraint_array[c];
    }

    inline void add_constraint(ConstraintBase* _constr) {
        constraint_array.push_back(std::shared_ptr<ConstraintBase>(_constr));
    }

    inline void add_constraint(std::shared_ptr<ConstraintBase>&& _constr) {
        constraint_array.push_back(std::shared_ptr<ConstraintBase>(_constr));
    }

    inline void delete_constraint(int c) {
        if (constraint_array[c]->get_can_be_deleted()) {
            constraint_array[c].reset();
        }
    }

    inline void delete_constraints(int first, int nb_del) {
        auto it = constraint_array.begin() + first;
        constraint_array.erase(it, it + nb_del);
    }
};

class BddCoeff : public VariableKeyBase {
   private:
    int    row;
    double coeff;
    double value;

   public:
    BddCoeff(int    _j,
             int    _t,
             double _coeff,
             double _value = 0.0,
             int    _row = -1,
             bool   _high = true,
             bool   _root = false)
        : VariableKeyBase(_j, _t, _high, _root),
          row(_row),
          coeff(_coeff),
          value(_value){};
    // BddCoeff() = default;
    ~BddCoeff() override = default;
    BddCoeff(const BddCoeff&) = default;
    BddCoeff& operator=(const BddCoeff&) = default;
    BddCoeff(BddCoeff&& op) = default;
    BddCoeff& operator=(BddCoeff&& op) = default;

    inline double get_coeff() { return coeff; }

    [[nodiscard]] inline double get_value() const { return value; }

    inline void set_value(double _value) { value = _value; }

    inline void set_row(int _row) { row = _row; }

    inline int get_row() { return row; }

    friend bool operator==(const BddCoeff& lhs, const BddCoeff& rhs) {
        return lhs.get_j() == rhs.get_j() && lhs.get_t() == rhs.get_t() &&
               lhs.get_high() == rhs.get_high();
    };

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
}  // namespace std

class GenericData {
   private:
    using coeff_hash_table = std::unordered_map<BddCoeff, double>;
    coeff_hash_table coeff;

    static constexpr double EPS_GENERIC_DATA = 1e-6;

   public:
    GenericData() = default;
    ~GenericData() = default;
    GenericData(GenericData&&) = default;  // movable and noncopyable
    GenericData(const GenericData&) = default;
    GenericData& operator=(GenericData&&) = default;
    GenericData& operator=(const GenericData&) = default;

    coeff_hash_table::iterator find(const BddCoeff& key) {
        return coeff.find(key);
    }

    coeff_hash_table::iterator end() { return coeff.end(); }

    void add_coeff_hash_table(int _j, int _t, bool _high, double _coeff) {
        BddCoeff key(_j, _t, _coeff, 0.0, -1, _high);

        auto it = coeff.find(key);
        if (it == coeff.end()) {
            coeff[key] = _coeff;
        } else {
            coeff[key] += _coeff;
        }
    }

    void list_coeff() {
        for (auto& it : coeff) {
            fmt::print("{} ({},{},{})", it.second, it.first.get_j(),
                       it.first.get_t(), it.first.get_high());
        }
        fmt::print("\n");
    }

    friend bool operator==(const GenericData& lhs, const GenericData& rhs) {
        if (lhs.coeff.size() != rhs.coeff.size()) {
            return false;
        }

        for (auto& it1 : lhs.coeff) {
            auto it2 = rhs.coeff.find(it1.first);
            if (it2 == rhs.coeff.end()) {
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

    double get_var_coeff(VariableKeyBase* key) override {
        auto* aux = static_cast<BddCoeff*>(key);
        auto  it = data->find(*aux);
        if (it == data->end()) {
            return 0.0;
        } else {
            return (*it).second;
        }
    }

    friend bool operator!=(const ConstraintGeneric& lhs,
                           const ConstraintGeneric& rhs) {
        auto& lhs_ptr = lhs.data;
        auto& rhs_ptr = rhs.data;
        return !(*lhs_ptr == *rhs_ptr);
    }

    friend bool operator==(const ConstraintGeneric& lhs,
                           const ConstraintGeneric& rhs) {
        auto& lhs_ptr = lhs.data;
        auto& rhs_ptr = rhs.data;
        return (*lhs_ptr == *rhs_ptr);
    }

    void list_coeff() { data->list_coeff(); }
};

template <typename T = BddCoeff>
class OriginalConstraint {
   private:
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

    inline std::list<std::shared_ptr<T>>* get_coeff_list() {
        return &coeff_list;
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

    inline void set_constraint(const std::shared_ptr<ConstraintBase>& _constr) {
        constr = _constr;
    }

    void clear_coeff() { coeff_list.clear(); }
};

template <typename T = BddCoeff>
class OriginalModel {
   private:
    std::vector<OriginalConstraint<T>> constraint_array;

   public:
    explicit OriginalModel(const ReformulationModel& model)
        : constraint_array(model.get_nb_constraints()) {
        const int nb_constraints = model.get_nb_constraints();

        for (int i = 0; i < nb_constraints; i++) {
            constraint_array[i].set_constraint(model.get_constraint_ptr(i));
        }
    }

    ~OriginalModel() = default;
    OriginalModel(OriginalModel&& op) noexcept =
        default;  // movable and noncopyable
    OriginalModel& operator=(OriginalModel&& op) noexcept = default;
    OriginalModel<T>(const OriginalModel<T>&) = default;
    OriginalModel<T>& operator=(const OriginalModel<T>&) = default;

    void add_coeff_list(int c, std::shared_ptr<T> coeff) {
        constraint_array[c].add_coeff_to_list(coeff);
    }

    ConstraintBase* get_constraint(int c) const {
        return constraint_array[c].get_constr();
    }

    inline std::list<std::shared_ptr<T>>* get_coeff_list(int c) const {
        return constraint_array[c].get_coeff_list();
    }

    inline void add_constraint(const std::shared_ptr<ConstraintBase>& _constr) {
        constraint_array.push_back(OriginalConstraint<>(_constr));
    }

    inline size_t get_nb_constraints() { return constraint_array.size(); }

    inline void delete_constraints(int first, int nb_del) {
        auto it = constraint_array.begin() + first;
        constraint_array.erase(it, it + nb_del);
    }

    void clear_all_coeff() {
        for (auto& it : constraint_array) {
            it.clear_coeff();
        }
    }
};

#endif