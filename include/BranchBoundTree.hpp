#ifndef __BRANCHBOUNDTREE_H__
#define __BRANCHBOUNDTREE_H__
#include <memory>
#include "branch-and-bound/btree.h"
#include "wctparms.h"
#include "wctprivate.h"

class BranchBoundTree {
   public:
    explicit BranchBoundTree(NodeData* data,
                             int       probType = 0,
                             bool      isIntProb = true);
    BranchBoundTree(BranchBoundTree&&) = default;
    BranchBoundTree& operator=(BranchBoundTree&&) = default;
    BranchBoundTree& operator=(const BranchBoundTree&) = delete;
    BranchBoundTree(const BranchBoundTree&) = delete;
    ~BranchBoundTree() = default;
    BTree* get_ptr_tree() { return tree.get(); }

   private:
    std::unique_ptr<BTree> tree;
};

#endif  // __BRANCHBOUNDTREE_H__