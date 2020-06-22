#ifndef _MODEL_INTERFACE
#define _MODEL_INTERFACE


#include "wctparms.h"
#include <memory>
#include <vector>
#include <list>
#include <NodeId.hpp>
class VariableKeyBase {
    private:
        int j;
        int t;
        bool root;


    public:
        VariableKeyBase(int _j, int _t, bool _root = false) : j(_j), t(_t), root(_root) {

        }

        VariableKeyBase() : j(-1), t(-1), root(false) {

        }
        inline int get_j() {
            return j;
        }

        inline int get_t() {
            return t;
        }

        inline bool get_root() {
            return root;
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
            if (key->get_j() == job) {
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
        if(key->get_t() == 0 && key->get_root()) {
            return -1.0;
        }

        return 0.0;
    }

};

class ReformulationModel {
private:
    std::vector<std::shared_ptr<ConstraintBase>> constraint_array;
    int nb_constraints;
public:
    ReformulationModel(int nb_assignments, int nb_machines);
    ~ReformulationModel();
    ReformulationModel(ReformulationModel && op) noexcept;              // movable and noncopyable
    ReformulationModel& operator=(ReformulationModel && op) noexcept;

    inline int get_nb_constraints() const {
        return nb_constraints;
    };

    inline ConstraintBase* get_constraint(int c) const {
        return constraint_array[c].get();
    };

};

class BddCoeff : public VariableKeyBase {
    private:
    int row;
    double coeff;
    bool high;

    public:
    BddCoeff(int _j, int _t, double _coeff, int _row = -1, bool _high = true, bool _root = false) : 
        VariableKeyBase(_j, _t, _root),
        row(_row),
        coeff(_coeff),
        high(_high) { } ;
    ~BddCoeff() = default;
    BddCoeff(BddCoeff && op) = default;              // movable and noncopyable
    BddCoeff& operator=(BddCoeff && op) = default;

    inline double get_coeff() {
        return coeff;
    } 

    inline bool get_high() {
        return high;
    }

    inline int get_row() {
        return row;
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


};


#endif