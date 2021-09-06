#ifndef __BRANCHBOUNDTREE_H__
#define __BRANCHBOUNDTREE_H__
#include <memory>                        // for unique_ptr
#include "branch-and-bound/TreeStats.h"  // for TreeStats
#include "branch-and-bound/btree.h"      // for BTree

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

    [[nodiscard]] double get_UB() const { return tree->getGlobalUB(); }
    [[nodiscard]] double get_LB() const { return tree->getGlobalLB(); }
    [[nodiscard]] int    get_nb_nodes_explored() const {
        return tree->tStats->get_states_explored();
    }

   private:
    std::unique_ptr<BTree> tree;
};

#endif  // __BRANCHBOUNDTREE_H__