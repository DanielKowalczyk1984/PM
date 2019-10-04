#ifndef NODE_BRANCH_ID_HPP
#define NODE_BRANCH_ID_HPP

struct NodeBranchId {
    size_t col;
    int row;
    int val;

    NodeBranchId() {
    }

    NodeBranchId(int row, size_t col, int val) :
            col(col), row(row), val(val) {
    }
};

#endif // NODE_BRANCH_ID_HPP
