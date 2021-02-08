#include "BranchBoundTree.hpp"
#include <fmt/core.h>
#include <boost/timer/timer.hpp>
#include <memory>
#include "BranchNode.hpp"
#include "branch-and-bound/bfstree.h"
#include "branch-and-bound/brfstree.h"
#include "branch-and-bound/cbfstree.h"
#include "branch-and-bound/dfstree.h"
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

    auto node = std::make_unique<BranchNodeBase>(root, true);
    tree->set_global_ub(double(root->opt_sol->tw));
    tree->set_retain_states(false);
    tree->set_time_limit(parms->branchandbound);
    node->set_id(root->id);
    node->set_depth(root->depth);
    fmt::print("{0:^10}|{1:^30}|{2:^30}|{3:^10}|{4:^10}|\n", "Nodes",
               "Current Node", "Objective Bounds", "Branch", "Work");
    fmt::print(
        "{0:^5}{1:^5}|{2:^10}{3:^10}{10:^10}|{4:>10}{5:>10}{6:>10}|{7:>5}{8:>5}"
        "|{9:>5}{11:>5}\n",
        "Expl", "Unex", "Obj", "Depth", "Primal", "Dual", "Gap", "Job", "Time",
        "Time", "Size", "Iter");

    tree->process_state(std::move(node), true);
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
