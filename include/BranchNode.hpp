#ifndef __BRANCHNODE_H__
#define __BRANCHNODE_H__

#include <limits>
#include <memory>
#include "branch-and-bound/state.h"
#include "wctprivate.h"

class BranchNodeBase : public State {
    struct BranchCand {
        double score{EPS_BRANCH};
        int    job{-1};

        BranchCand() = default;

        bool operator<(const BranchCand& other) const {
            return (score < other.score);
        };
    };

   public:
    explicit BranchNodeBase(NodeData* pd, bool isRoot = false);
    BranchNodeBase(BranchNodeBase&&) = default;
    BranchNodeBase(const BranchNodeBase&) = default;
    BranchNodeBase& operator=(BranchNodeBase&&) = default;
    BranchNodeBase& operator=(const BranchNodeBase&) = default;
    ~BranchNodeBase() override { nodedata_free(pd); };
    static constexpr double EPS_BRANCH = 1e-4;
    static constexpr double ERROR = 1e-12;

    void branch(BTree* bt) final;
    void compute_bounds(BTree* bt) final;
    void assess_dominance(State* otherState) final;
    bool is_terminal_state() final;
    void apply_final_pruning_tests(BTree* bt) final;
    void update_data(double upper_bound) final;
    // std::unique_ptr<State> clone() { return nullptr; };  // "copy
    // constructor"
    void print() const final{};
    bool operator<(const State& other) final { return false; };

   private:
    NodeData* pd;

    static constexpr double TargetBrTimeValue = 0.5;
    static constexpr int    NumStrBrCandidates = 16;
    static constexpr double IntegerTolerance = 1e-3;
};

#endif  // __BRANCHNODE_H__