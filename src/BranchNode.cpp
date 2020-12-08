#include "BranchNode.hpp"
#include "solver.h"
#include "wct.h"
#include "wctprivate.h"

BranchNodeBase::BranchNodeBase(NodeData* _pd, bool _isRoot)
    : State(_isRoot),
      pd(_pd) {}

void BranchNodeBase::branch(BTree* bt) {}

void BranchNodeBase::computeBounds(BTree* bt) {
    // build_rmp(pd);
    // solve_relaxation(pd);
    compute_lower_bound(pd);
    lowerBound = pd->LP_lower_bound;
    objValue = pd->opt_sol->tw;
}

void BranchNodeBase::assessDominance(State* otherState) {}

bool BranchNodeBase::isTerminalState() {
    construct_lp_sol_from_rmp(pd);
    return false;
}

void BranchNodeBase::applyFinalPruningTests(BTree* bt) {}

extern "C" {

BranchNodeBase* new_branch_node(int _isRoot, NodeData* data) {
    return new BranchNodeBase(data, _isRoot);
};

void delete_branch_node(BranchNodeBase* node) {
    if (node) {
        delete node;
    }
}

size_t call_getDepth(BranchNodeBase* state) {
    return state->getDepth();
};
int call_getDomClassID(BranchNodeBase* state) {
    return state->getDomClassID();
};
double call_getObjValue(BranchNodeBase* state) {
    return state->getObjValue();
};
double call_getLB(BranchNodeBase* state) {
    return state->getLB();
};
double call_getUB(BranchNodeBase* state) {
    return state->getUB();
};
int call_getID(BranchNodeBase* state) {
    return state->getID();
}
int call_getParentID(BranchNodeBase* state) {
    return state->getParentID();
}
void call_setID(BranchNodeBase* state, int i) {
    state->setID(i);
}
bool call_isDominated(BranchNodeBase* state) {
    return state->isDominated();
};
bool call_wasProcessed(BranchNodeBase* state) {
    return state->wasProcessed();
};
}
