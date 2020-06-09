#include "ModelInterface.hpp"
#include <memory>

ReformulationModel::ReformulationModel(int nb_assignments, int nb_machines) : constraint_array(nb_assignments + 1, nullptr) {
    int counter = 0;
    for(auto &it: constraint_array) {
        it=  std::make_shared<ConstraintAssignment>(counter++);
    }

    constraint_array[counter] = std::make_shared<ConstraintConvex>((double) nb_machines);

}

ReformulationModel::~ReformulationModel() = default;
ReformulationModel::ReformulationModel(ReformulationModel &&) noexcept = default;
ReformulationModel& ReformulationModel::operator=(ReformulationModel &&) noexcept = default;