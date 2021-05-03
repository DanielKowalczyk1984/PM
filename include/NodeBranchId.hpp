#ifndef NODE_BRANCH_ID_HPP
#define NODE_BRANCH_ID_HPP
#include <cstddef>

struct NodeBranchId {
    size_t col;
    size_t row;
    int    val;

    NodeBranchId() = default;

    NodeBranchId(size_t _row, size_t _col, int _val)
        : col(_col),
          row(_row),
          val(_val) {}
};

#endif  // NODE_BRANCH_ID_HPP
