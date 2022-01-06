// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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