#include "BranchHistory.hpp"
#include <range/v3/view/zip.hpp>

double BranchHistory::update_pseudocost(
    const std::function<double(double, double)>& function,
    std::array<double, 2>&                       gain) {
    for (auto&& [s, e, p, g] : ranges::views::zip(sigma, eta, psi, gain)) {
        s += g;
        e++;
        p = s / e;
    }

    score = function(gain[0], gain[1]);

    return score;
}
std::array<double, 2> BranchHistory::psi_average = {1.0, 1.0};
std::size_t           BranchHistory::nb_histroy = 0UL;