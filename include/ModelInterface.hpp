#ifndef _MODEL_INTERFACE
#define _MODEL_INTERFACE


#include "wctparms.h"
#include <memory>
#include <vector>
class VariableKeyBase {
    private:
        int j;
        int t;

    public:
        VariableKeyBase(int _j, int _t) : j(_j), t(_t) {

        }

        VariableKeyBase() : j(-1), t(-1) {

        }
        inline int get_j() {
            return j;
        }

        inline int get_t() {
            return t;
        }

        inline void set_j(int _j) {
            j = _j;
        }

        inline void set_t(int _t) {
            t = _t;
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
        return -1.0;
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

    inline int get_nb_constraints() {
        return nb_constraints;
    };

    inline ConstraintBase* get_constraint(int c) {
        return constraint_array[c].get();
    };

};


#endif