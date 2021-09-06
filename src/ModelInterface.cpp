#include "ModelInterface.hpp"
#include <memory>                                // for make_shared, shared_ptr
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=, operator-
#include <range/v3/iterator/diffmax_t.hpp>       // for operator<=
#include <range/v3/range/conversion.hpp>         // for operator|, to, to_co...
#include <range/v3/view/iota.hpp>                // for iota_view, iota, iot...
#include <range/v3/view/transform.hpp>           // for transform_view, tran...
#include <range/v3/view/view.hpp>                // for operator|, view_closure
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
