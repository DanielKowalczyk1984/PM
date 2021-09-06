#include "BranchBoundTree.hpp"
#include <memory>                       // for make_unique, unique_ptr
#include <utility>                      // for move
#include "BranchNode.hpp"               // for BranchNodeBase
#include "NodeData.h"                   // for NodeData
#include "Parms.h"                      // for Parms, bb_bfs_strategy, bb_br...
#include "Solution.hpp"                 // for Sol
#include "branch-and-bound/bfstree.h"   // for BFSTree
#include "branch-and-bound/brfstree.h"  // for BrFSTree
#include "branch-and-bound/cbfstree.h"  // for CBFSTree
#include "branch-and-bound/dfstree.h"   // for DFSTree
#include "branch-and-bound/state.h"     // for State

BranchBoundTree::BranchBoundTree(std::unique_ptr<NodeData> root,
                                 int                       _probType,
                                 bool                      _isIntProb) {
    auto& parms = root->parms;
    switch (parms.bb_explore_strategy) {
        case min_bb_explore_strategy:
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

    tree->set_global_ub(double(root->opt_sol.tw));
    tree->set_retain_states(false);
    tree->set_final_test_usage(parms.pruning_test);
    tree->set_time_limit(static_cast<int>(parms.branching_cpu_limit));
    tree->set_node_limit(parms.bb_node_limit);
    tree->set_only_root_node(parms.bb_node_limit == 1);
    auto aux = root->depth;
    auto node = std::make_unique<BranchNodeBase>(std::move(root), true);
    node->set_depth(aux);
    tree->process_state(std::move(node), true);
}
