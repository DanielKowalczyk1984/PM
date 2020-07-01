#ifndef _MODEL_INTERFACE
#define _MODEL_INTERFACE


#include <bits/c++config.h>
#include "wctparms.h"
#include <boost/container_hash/hash_fwd.hpp>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include <list>
#include <NodeId.hpp>
class VariableKeyBase {
    private:
        int j;
        int t;
        bool high;
        bool root;


    public:
        VariableKeyBase(int _j, int _t, bool _high = true, bool _root = false) : j(_j), t(_t), high(_high),  root(_root) {

        }

        VariableKeyBase() : j(-1), t(-1), high(false), root(false) {

        }
        inline int get_j() const {
            return j;
        }

        inline int get_t() const {
            return t;
        }

        inline bool get_root() {
            return root;
        }

        inline bool get_high() const {
            return high;
        }

        inline void set_j(int _j) {
            j = _j;
        }

        inline void set_t(int _t) {
            t = _t;
        }

        inline void set_root(bool _root) {
            root = _root;
        }

        inline void set_high(bool _high) {
            high = _high;
        }

        VariableKeyBase(const VariableKeyBase&) = default;
        VariableKeyBase(VariableKeyBase&&) = default;
        VariableKeyBase& operator=(const VariableKeyBase&) = default;
        VariableKeyBase& operator=(VariableKeyBase&&) = default;
        virtual ~VariableKeyBase() = default;

};

class ConstraintBase {
    protected:
        char sense;
        double rhs;
        bool can_be_deleted;
    
    public:
        inline double get_rhs() {
            return rhs;
        }

        inline char get_sense() {
            return sense;
        }

        inline bool get_can_be_deleted () {
            return can_be_deleted;
        }

        ConstraintBase(const ConstraintBase&) = default;
        ConstraintBase(ConstraintBase&&) = default;
        ConstraintBase& operator=(const ConstraintBase&) = default;
        ConstraintBase& operator=(ConstraintBase&&) = default;
        virtual ~ConstraintBase() = default;

        ConstraintBase(char _sense, double _rhs) : sense(_sense), rhs(_rhs), can_be_deleted(false) {

        }

        virtual double get_var_coeff(VariableKeyBase*) = 0;
};

class ConstraintAssignment : public ConstraintBase {
    private:
        int job;
    public:
        ConstraintAssignment(int _job) : ConstraintBase('>', 1.0), job(_job) {

        }

        double get_var_coeff(VariableKeyBase *key) {
            if (key->get_j() == job && key->get_high()) {
                return 1.0;
            }

            return 0.0;
        }

        ConstraintAssignment(const ConstraintAssignment&) = default;
        ConstraintAssignment(ConstraintAssignment&&) = default;
        ConstraintAssignment& operator=(const ConstraintAssignment&) = default;
        ConstraintAssignment& operator=(ConstraintAssignment&&) = default;
        virtual ~ConstraintAssignment() = default;

};

class ConstraintConvex : public ConstraintBase {
    public:
    ConstraintConvex(double _rhs): ConstraintBase('>',_rhs) {

    }

    double get_var_coeff(VariableKeyBase *key) {
        if(key->get_t() == 0 ) {
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
    ~ReformulationModel();
    ReformulationModel(ReformulationModel && op) noexcept;              // movable and noncopyable
    ReformulationModel& operator=(ReformulationModel && op) noexcept;

    inline int get_nb_constraints() const {
        return constraint_array.size();
    };

    inline ConstraintBase* get_constraint(int c) const {
        return constraint_array[c].get();
    };

    inline void add_constraint(ConstraintBase* _constr) {
        constraint_array.push_back(std::shared_ptr<ConstraintBase>(_constr));
    } 

};

class BddCoeff : public VariableKeyBase {
    private:
    int row;
    double coeff;
    double value;

    public:
    BddCoeff(int _j, int _t, double _coeff, double _value = 0.0, int _row = -1, bool _high = true, bool _root = false) : 
        VariableKeyBase(_j, _t, _high, _root),
        row(_row),
        coeff(_coeff),
        value(_value) { } ;
    BddCoeff() = default;
    ~BddCoeff() = default;
    BddCoeff(const BddCoeff&) = default;
    BddCoeff& operator=(const BddCoeff&) = default;
    BddCoeff(BddCoeff && op) = default;
    BddCoeff& operator=(BddCoeff && op) = default;

    inline double get_coeff() {
        return coeff;
    }

    inline double get_value() {
        return value;
    }

    inline void set_value(double _value) {
        value = _value;
    }

    inline int get_row() {
        return row;
    }

    friend bool operator==(const BddCoeff& lhs, const BddCoeff & rhs) {
        return lhs.get_j() == rhs.get_j() && lhs.get_t() == rhs.get_t() && lhs.get_high() == rhs.get_high();
    }
};

namespace std {
template <>
struct hash<BddCoeff> {
    std::size_t operator()(BddCoeff const& s) const noexcept {
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
    
    typedef std::unordered_map<BddCoeff, double> coeff_hash_table;
    coeff_hash_table coeff;


    public:

    GenericData() {} 
    ~GenericData() = default;
    GenericData(GenericData && op) = default;              // movable and noncopyable
    GenericData& operator=(GenericData && op) = default;

    coeff_hash_table::iterator find(const BddCoeff & key) {
        return coeff.find(key);
    }

    coeff_hash_table::iterator end() {
        return coeff.end();
    }

    void add_coeff_hash_table(int _j, int _t, bool _high, double _coeff) {
        BddCoeff key(_j,_t, _coeff, 0.0, -1, _high);

        auto it = coeff.find(key);
        if (it == coeff.end()) {
            coeff[key] = _coeff;
        } else {
            coeff[key] += _coeff;
        }
    }

    void list_coeff() {
        for(auto &it: coeff) {
            std::cout << it.second << " ";
        }
        std::cout << "\n";
    }

    friend bool operator==(const GenericData& lhs, const GenericData& rhs) {
        if (lhs.coeff.size() != rhs.coeff.size()) {
            return false;
        }

        for(auto& it1: lhs.coeff) {
            auto it2 = rhs.coeff.find(it1.first);
            if(it2 == rhs.coeff.end()) {
                return false;
            }

            if(fabs(it1.second - (*it2).second) > 1e-6) {
                return false;
            }
        }

        return true;
    };

    friend bool operator!=(const GenericData& lhs, GenericData & rhs) { return !(lhs == rhs); }
};
class ConstraintGeneric : public ConstraintBase {
    private:
    std::unique_ptr<GenericData> data;
    
    public:
    ConstraintGeneric(GenericData* _data, double _rhs, char _sense = '>') : ConstraintBase(_sense, _rhs), data(_data) {
    
    }

    ConstraintGeneric(double _rhs, char _sense = '>') : ConstraintBase(_sense, _rhs) {

    }

    ~ConstraintGeneric() = default;
    ConstraintGeneric(ConstraintGeneric && op) = default;              // movable and noncopyable
    ConstraintGeneric& operator=(ConstraintGeneric && op) = default;

    double get_var_coeff(VariableKeyBase *key) override {
        BddCoeff* aux = static_cast<BddCoeff*>(key);
        auto it = data->find(*aux);
        if (it == data->end()) {
            return 0.0;
        } else {
            return (*it).second;
        }
    }

    void list_coeff() {
        data->list_coeff();
    }

};


template<typename T = BddCoeff>
class OriginalConstraint {
    private:
    ConstraintBase* constr;
    std::list<std::shared_ptr<T>> coeff_list;

    public:
    OriginalConstraint(ConstraintBase* _constr): constr(_constr) { } ;
    OriginalConstraint(): constr(nullptr) { } ;
    ~OriginalConstraint() = default;
    OriginalConstraint(OriginalConstraint && op) = default;              // movable and noncopyable
    OriginalConstraint& operator=(OriginalConstraint && op) = default;

    inline std::list<std::shared_ptr<T>>* get_coeff_list() {
        return &coeff_list;
    }

    inline ConstraintBase* get_constr() {
        return constr;
    }

    inline void add_coeff_to_list(std::shared_ptr<T> _coeff) {
        coeff_list.push_back(_coeff);
    }

    inline void set_constraint(ConstraintBase* _constr) {
        constr = _constr;
    }
};

template<typename T = BddCoeff>
class OriginalModel{
    private:
    std::vector<OriginalConstraint<T>> constraint_array;

    public:
    OriginalModel(const ReformulationModel& model) : constraint_array(model.get_nb_constraints()) {
        const int nb_constraints = model.get_nb_constraints();

        for(int i = 0; i < nb_constraints; i++) {
            constraint_array[i].set_constraint(model.get_constraint(i));
        }
    }

    ~OriginalModel() = default;
    OriginalModel(OriginalModel && op) = default;              // movable and noncopyable
    OriginalModel& operator=(OriginalModel && op) = default;
    
    void add_coeff_list(int c,  std::shared_ptr<T> coeff) {
        constraint_array[c].add_coeff_to_list(coeff);
    }

    ConstraintBase* get_constraint(int c) {
        return constraint_array[c].get_constr();
    }

    inline std::list<std::shared_ptr<T>>* get_coeff_list(int c) {
        return constraint_array[c].get_coeff_list();
    }

    inline void add_constraint(ConstraintBase* _constr) {
        constraint_array.push_back(OriginalConstraint<>(_constr));
    }

    inline size_t get_nb_constraints() {
        return constraint_array.size();
    }


};


#endif