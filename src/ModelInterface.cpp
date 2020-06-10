#include "ModelInterface.hpp"
#include <memory>

ReformulationModel::ReformulationModel(int nb_assignments, int nb_machines) : constraint_array(nb_assignments + 1, nullptr) {
    nb_constraints = 0;
    for(int i = 0; i < nb_assignments; i++) {
        constraint_array[i] = std::make_shared<ConstraintAssignment>(i);
        nb_constraints++;
    }

    constraint_array[nb_constraints] = std::make_shared<ConstraintConvex>((double) -nb_machines);

}

ReformulationModel::~ReformulationModel() = default;
ReformulationModel::ReformulationModel(ReformulationModel &&) noexcept = default;
ReformulationModel& ReformulationModel::operator=(ReformulationModel &&) noexcept = default;