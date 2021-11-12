// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
