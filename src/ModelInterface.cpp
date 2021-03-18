#include "ModelInterface.hpp"
#include <memory>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

ReformulationModel::ReformulationModel(int nb_assignments, int nb_machines)
    : std::vector<std::shared_ptr<ConstraintBase>>(
          ranges::views::iota(0, nb_assignments) |
          ranges::views::transform(
              [](int i) { return std::make_shared<ConstraintAssignment>(i); }) |
          ranges::to<std::vector<std::shared_ptr<ConstraintBase>>>()) {
    (*this).push_back(
        std::make_shared<ConstraintConvex>(static_cast<double>(-nb_machines)));
}
