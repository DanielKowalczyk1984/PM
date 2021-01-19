#ifndef __BRANCHNODE_H__
#define __BRANCHNODE_H__

#include <limits>
#include "branch-and-bound/state.h"
#include "wctprivate.h"

class BranchNodeBase : public State {
    struct BranchCand {
        double score{EPS};
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
    static constexpr double EPS = 1e-4;
    static constexpr double ERROR = 1e-12;

    void   branch(BTree* bt) final;
    void   computeBounds(BTree* bt) final;
    void   assessDominance(State* otherState) final;
    bool   isTerminalState() final;
    void   applyFinalPruningTests(BTree* bt) final;
    State* clone() final { return nullptr; };  // "copy constructor"
    void   print() const final{};
    bool   operator<(const State& other) final { return false; };

   private:
    NodeData* pd;

    static constexpr double TargetBrTimeValue = 0.5;
    static constexpr int    NumStrBrCandidates = 16;
    static constexpr double IntegerTolerance = 1e-3;
};

#endif  // __BRANCHNODE_H__