#ifndef NODE_BRANCH_ID_HPP
#define NODE_BRANCH_ID_HPP
#include <stddef.h>

struct NodeBranchId {
    size_t col;
    int    row;
    int    val;

    NodeBranchId() {}

    NodeBranchId(int _row, size_t _col, int _val)
        : col(_col),
          row(_row),
          val(_val) {}
};

#endif  // NODE_BRANCH_ID_HPP
