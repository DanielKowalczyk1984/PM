#ifndef NODE_BRANCH_ID_HPP
#define NODE_BRANCH_ID_HPP
#include <cstddef>

struct NodeBranchId {
    size_t col{};
    size_t row{};
    size_t val{};

    NodeBranchId() = default;

    NodeBranchId(size_t _row, size_t _col, size_t _val)
        : col(_col),
          row(_row),
          val(_val) {}
} __attribute__((packed)) __attribute__((aligned(32)));

#endif  // NODE_BRANCH_ID_HPP
