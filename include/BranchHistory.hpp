#ifndef __BRANCHHISTORY_H__
#define __BRANCHHISTORY_H__

#include <array>
#include <functional>
#include "NodeId.hpp"
class BranchHistory {
   public:
    BranchHistory();
    BranchHistory(BranchHistory&&) = default;
    BranchHistory(const BranchHistory&) = default;
    BranchHistory& operator=(BranchHistory&&) = default;
    BranchHistory& operator=(const BranchHistory&) = default;
    ~BranchHistory() = default;

    double update_pseudocost(const std::function<double(double, double)>&,
                             std::array<double, 2>& gain);

    double initial_score(const std::function<double(double, double)>&);

    void compute_score(const std::function<double(double, double)>&,
                       double left,
                       double right);

   private:
    std::array<double, 2> sigma{};
    std::array<double, 2> psi{};
    std::array<int, 2>    eta{};
    double                score{};

    static std::array<double, 2> psi_average;
    static std::size_t           nb_histroy;
};

#endif  // __BRANCHHISTORY_H__