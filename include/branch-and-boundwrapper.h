#ifndef __BRANCH_AND_BOUNDWRAPPER_H__
#define __BRANCH_AND_BOUNDWRAPPER_H__

#include "wctprivate.h"

#ifdef __cplusplus
extern "C" {

#endif  // __cplusplus

typedef struct BranchNodeBase  BranchNode;
typedef struct BranchBoundTree BranchBoundTree;

BranchNode* new_branch_node(int _isRoot, NodeData* data);
void        delete_branch_node(BranchNode* node);
size_t      call_getDepth(BranchNode* state);
int         call_getDomClassID(BranchNode* state);
double      call_getObjValue(BranchNode* state);
double      call_getLB(BranchNode* state);
double      call_getUB(BranchNode* state);
int         call_getID(BranchNode* state);
int         call_getParentID(BranchNode* state);
void        call_setID(BranchNode* state, int i);
int         isDominated(BranchNode* state);
int         wasProcessed(BranchNode* state);

BranchBoundTree* new_branch_bound_tree(NodeData* data,
                                       int       _probtype,
                                       int       _isIntProb);
void             delete_branch_bound_tree(BranchBoundTree* tree);

#ifdef __cplusplus
}
#endif  // __cplusplus
#endif  // __BRANCH-AND-BOUNDWRAPPER_H__
