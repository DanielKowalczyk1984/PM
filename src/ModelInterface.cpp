#include "ModelInterface.hpp"
#include <memory>
#include <utility>

ReformulationModel::ReformulationModel(int nb_assignments, int nb_machines)
    : constraint_array(nb_assignments + 1, nullptr) {
    for (int i = 0; i < nb_assignments; i++) {
        constraint_array[i] = std::make_shared<ConstraintAssignment>(i);
    }

    constraint_array[nb_assignments] =
        std::make_shared<ConstraintConvex>(static_cast<double>(-nb_machines));

    // GenericData* data = new GenericData();

    // std::vector<std::pair<int, int>> coeff_data = { {3, 71},{3,110},{4,71},
    // {4,152} };

    // for(auto &it: coeff_data) {
    //     (*data)[BddCoeff(it.first,it.second,0.0,0.0)] = -1.0;
    // }

    // std::shared_ptr<ConstraintGeneric> constr(new
    // ConstraintGeneric(data,-1.0,'>'));

    // constraint_array.push_back(constr);
}

ReformulationModel::~ReformulationModel() = default;
ReformulationModel::ReformulationModel(ReformulationModel&&) noexcept = default;
ReformulationModel& ReformulationModel::operator=(
    ReformulationModel&&) noexcept = default;
