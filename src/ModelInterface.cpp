#include "ModelInterface.hpp"
#include <memory>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

ReformulationModel::ReformulationModel(size_t nb_assignments,
                                       size_t nb_machines)
    : std::vector<std::shared_ptr<ConstraintBase>>(
          ranges::views::iota(0UL, nb_assignments) |
          ranges::views::transform([](auto i) {
              return std::make_shared<ConstraintAssignment>(i);
          }) |
          ranges::to<std::vector<std::shared_ptr<ConstraintBase>>>()) {
    auto m = static_cast<double>(nb_machines);
    (*this).push_back(std::make_shared<ConstraintConvex>(-m));
}
