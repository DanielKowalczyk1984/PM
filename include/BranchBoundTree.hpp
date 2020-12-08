#ifndef __BRANCHBOUNDTREE_H__
#define __BRANCHBOUNDTREE_H__
#include <memory>
#include "branch-and-bound/btree.h"
#include "wctparms.h"
#include "wctprivate.h"

class BranchBoundTree {
   public:
    BranchBoundTree(NodeData* data, int probType = 0, bool isIntProb = 1);
    BranchBoundTree(BranchBoundTree&&) = default;
    BranchBoundTree& operator=(BranchBoundTree&&) = default;
    ~BranchBoundTree(){};

   private:
    std::unique_ptr<BTree> tree;
};

#endif  // __BRANCHBOUNDTREE_H__