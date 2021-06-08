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
    void print(const BTree* bt) const override;

    [[nodiscard]] NodeData* get_data_ptr() const { return pd.get(); }

   private:
    std::unique_ptr<NodeData> pd;

    static constexpr double ERROR = 1e-12;
    static constexpr double IntegerTolerance = 1e-3;
};

struct BranchCand {
    double score{EPS_BRANCH};
    int    job{-1};
    int    t{-1};

    // std::unique_ptr<BranchNodeBase> left{};
    // std::unique_ptr<BranchNodeBase> right{};

    BranchCand() = default;

    BranchCand(double _score, int _job, int _t);

    static constexpr double EPS_BRANCH = 1e-4;
};
#endif  // __BRANCHNODE_H__