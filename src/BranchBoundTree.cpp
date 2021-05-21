#include "BranchBoundTree.hpp"
#include <memory>
#include "BranchNode.hpp"
#include "Parms.h"
#include "Statistics.h"
#include "branch-and-bound/bfstree.h"
#include "branch-and-bound/brfstree.h"
#include "branch-and-bound/cbfstree.h"
#include "branch-and-bound/dfstree.h"
#include "wctprivate.h"

BranchBoundTree::BranchBoundTree(std::unique_ptr<NodeData> root,
                                 int                       _probType,
                                 bool                      _isIntProb) {
    switch (root->parms.bb_explore_strategy) {
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
    tree->set_time_limit(root->parms.branching_cpu_limit);
    tree->set_node_limit(root->parms.bb_node_limit);
    auto aux = root->depth;
    auto node = std::make_unique<BranchNodeBase>(std::move(root), true);
    node->set_depth(aux);
    tree->process_state(std::move(node), true);
}
