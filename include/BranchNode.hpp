#ifndef __BRANCHNODE_H__
#define __BRANCHNODE_H__

#include "branch-and-bound/state.h"
#include "wctprivate.h"

class BranchNodeBase : public State {
   public:
    BranchNodeBase(NodeData* pd, bool isRoot = false);
    BranchNodeBase(BranchNodeBase&&) = default;
    BranchNodeBase(const BranchNodeBase&) = default;
    BranchNodeBase& operator=(BranchNodeBase&&) = default;
    BranchNodeBase& operator=(const BranchNodeBase&) = default;
    virtual ~BranchNodeBase(){};

    virtual void   branch(BTree* bt) final;
    virtual void   computeBounds(BTree* bt) final;
    virtual void   assessDominance(State* otherState) final;
    virtual bool   isTerminalState() final;
    virtual void   applyFinalPruningTests(BTree* bt) final;
    virtual State* clone() { return nullptr; };  // "copy constructor"
    virtual void   print() const {};
    virtual bool   operator<(const State& other) { return false; };

   private:
    NodeData* pd;
};

#endif  // __BRANCHNODE_H__