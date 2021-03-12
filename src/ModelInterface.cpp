#include "ModelInterface.hpp"
#include <memory>
#include <utility>

ReformulationModel::ReformulationModel(int nb_assignments, int nb_machines)
    : std::vector<std::shared_ptr<ConstraintBase>>(nb_assignments + 1,
                                                   nullptr) {
    for (int i = 0; i < nb_assignments; i++) {
        (*this)[i] = std::make_shared<ConstraintAssignment>(i);
    }

    (*this)[nb_assignments] =
        std::make_shared<ConstraintConvex>(static_cast<double>(-nb_machines));
}
