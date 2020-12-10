#include "BranchBoundTree.hpp"
#include <fmt/core.h>
#include <memory>
#include "BranchNode.hpp"
#include "branch-and-bound/bfstree.h"
#include "branch-and-bound/brfstree.h"
#include "branch-and-bound/cbfstree.h"
#include "branch-and-bound/dfstree.h"
#include "branch-and-boundwrapper.h"
#include "solution.h"
#include "wctparms.h"
#include "wctprivate.h"

BranchBoundTree::BranchBoundTree(NodeData* root,
                                 int       _probType,
                                 bool      _isIntProb) {
    Parms* parms = root->parms;
    switch (parms->bb_explore_strategy) {
        case min_bb_strategy:
            tree = std::make_unique<DFSTree>(_probType, _isIntProb);
            break;
        case bb_bfs_strategy:
            tree = std::make_unique<BFSTree>(_probType, _isIntProb);
            break;
        case bb_brfs_strategy:
            tree = std::make_unique<BrFSTree>(_probType, _isIntProb);
            break;
        case bb_cbfs_strategy:
            tree = std::make_unique<CBFSTree>(_probType, _isIntProb);
            break;
    }

    BranchNodeBase* node = new BranchNode(root, true);
    tree->setGlobalUB(double(root->opt_sol->tw));
    node->setID(root->id);
    node->setDepth(root->depth);
    tree->processState(node, true);
}

extern "C" {
BranchBoundTree* new_branch_bound_tree(NodeData* root,
                                       int       _probType,
                                       int       _isIntProb) {
    return new BranchBoundTree(root, _probType, _isIntProb);
}

void delete_branch_bound_tree(BranchBoundTree* tree) {
    if (tree) {
        delete tree;
    }
}

void call_branch_and_bound_explore(BranchBoundTree* tree) {
    BTree* ptr_tree = tree->get_ptr_tree();
    ptr_tree->explore();
}
}
