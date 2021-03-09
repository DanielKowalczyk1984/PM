#ifndef __BRANCHBOUNDTREE_H__
#define __BRANCHBOUNDTREE_H__
#include <memory>
#include "Parms.h"
#include "branch-and-bound/btree.h"

struct NodeData;
class BranchBoundTree {
   public:
    explicit BranchBoundTree(std::unique_ptr<NodeData> data,
                             int                       probType = 0,
                             bool                      isIntProb = true);
    BranchBoundTree(BranchBoundTree&&) = default;
    BranchBoundTree& operator=(BranchBoundTree&&) = default;
    BranchBoundTree& operator=(const BranchBoundTree&) = delete;
    BranchBoundTree(const BranchBoundTree&) = delete;
    ~BranchBoundTree() = default;
    void explore() { tree->explore(); }

   private:
    std::unique_ptr<BTree> tree;
};

#endif  // __BRANCHBOUNDTREE_H__