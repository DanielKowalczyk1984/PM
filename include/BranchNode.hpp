#ifndef __BRANCHNODE_H__
#define __BRANCHNODE_H__

#include <fmt/core.h>
#include <limits>
#include <memory>
#include "branch-and-bound/btree.h"
#include "branch-and-bound/state.h"
// #include "wctprivate.h"
class NodeData;
class BranchNodeBase : public State {
   public:
    explicit BranchNodeBase(std::unique_ptr<NodeData> pd, bool isRoot = false);
    BranchNodeBase(BranchNodeBase&&) = default;
    BranchNodeBase(const BranchNodeBase&) = delete;
    BranchNodeBase& operator=(BranchNodeBase&&) = default;
    BranchNodeBase& operator=(const BranchNodeBase&) = delete;
    ~BranchNodeBase() override = default;
    [[nodiscard]] std::unique_ptr<State> clone() const override {
        return nullptr;
    }

    void branch(BTree* bt) final;
    void compute_bounds(BTree* bt) final;
    void assess_dominance(State* otherState) final;
    bool is_terminal_state() final;
    void apply_final_pruning_tests(BTree* bt) final;
    void update_data(double upper_bound) final;
    // std::unique_ptr<State> clone() { return nullptr; };  // "copy
    // constructor"
    void print(const BTree* bt) const override;
    bool operator<(const State& other) final {
        return get_obj_value() < other.get_obj_value();
    };

   private:
    std::unique_ptr<NodeData> pd;

    static constexpr double ERROR = 1e-12;
    static constexpr double IntegerTolerance = 1e-3;
    static constexpr double TargetBrTimeValue = 0.3;
    static constexpr int    NumStrBrCandidates = 32;
};

struct BranchCand {
    double score{EPS_BRANCH};
    int    job{-1};
    int    t{-1};

    // std::unique_ptr<BranchNodeBase> left{nullptr};
    // std::unique_ptr<BranchNodeBase> right{nullptr};

    BranchCand() = default;

    bool operator<(const BranchCand& other) const {
        return (score < other.score);
    };
    static constexpr double EPS_BRANCH = 1e-4;
};

namespace std {
template <>
struct less<BranchCand> {
    constexpr bool operator()(auto const& lhs, auto const& rhs) {
        return rhs < lhs;  // or use boost::hash_combine
    }
};
}  // namespace std

#endif  // __BRANCHNODE_H__