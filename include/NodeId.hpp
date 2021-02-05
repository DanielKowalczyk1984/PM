#ifndef NODE_ID_HPP
#define NODE_ID_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <ostream>

int const NODE_ROW_BITS = 20;
int const NODE_ATTR_BITS = 1;
int const NODE_COL_BITS = 64 - NODE_ROW_BITS - NODE_ATTR_BITS;

int const NODE_ROW_OFFSET = NODE_COL_BITS + NODE_ATTR_BITS;
int const NODE_ATTR_OFFSET = NODE_COL_BITS;

uint64_t const NODE_ROW_MAX = (uint64_t(1) << NODE_ROW_BITS) - 1;
uint64_t const NODE_COL_MAX = (uint64_t(1) << NODE_COL_BITS) - 1;

uint64_t const NODE_ROW_MASK = NODE_ROW_MAX << NODE_ROW_OFFSET;
uint64_t const NODE_ATTR_MASK = uint64_t(1) << NODE_ATTR_OFFSET;

class NodeId {
    uint64_t code_;

    static constexpr size_t HASH_CONSTANT = 314159257;

   public:
    NodeId() = default;
    // {  // 'code_' is not initialized in the default constructor for
    //    // SPEED.
    // }

    NodeId(uint64_t code) : code_(code) {}

    NodeId(uint64_t row, uint64_t col)
        : code_((row << NODE_ROW_OFFSET) | col) {}

    NodeId(uint64_t row, uint64_t col, bool attr)
        : code_((row << NODE_ROW_OFFSET) | col) {
        setAttr(attr);
    }

    [[nodiscard]] size_t row() const { return code_ >> NODE_ROW_OFFSET; }

    [[nodiscard]] size_t col() const { return code_ & NODE_COL_MAX; }

    void setAttr(bool val) {
        if (val) {
            code_ |= NODE_ATTR_MASK;
        } else {
            code_ &= ~NODE_ATTR_MASK;
        }
    }

    [[nodiscard]] bool getAttr() const { return (code_ & NODE_ATTR_MASK) != 0; }

    [[nodiscard]] NodeId withoutAttr() const { return code_ & ~NODE_ATTR_MASK; }

    [[nodiscard]] bool hasEmpty() const { return code_ == 1 || getAttr(); }

    [[nodiscard]] uint64_t code() const { return code_ & ~NODE_ATTR_MASK; }

    [[nodiscard]] size_t hash() const { return code() * HASH_CONSTANT; }

    bool operator==(NodeId const& o) const { return code() == o.code(); }

    bool operator==(NodeId const&& o) const { return code() == o.code(); }

    bool operator!=(NodeId const& o) const { return !(*this == o); }

    bool operator<(NodeId const& o) const { return code() < o.code(); }

    bool operator>=(NodeId const& o) const { return !(*this < o); }

    bool operator>(NodeId const& o) const { return o < *this; }

    bool operator<=(NodeId const& o) const { return !(o < *this); }

    NodeId& operator=(const NodeId&) = default;
    NodeId(const NodeId&) = default;
    NodeId(NodeId&&) = default;
    NodeId& operator=(NodeId&&) = default;
    ~NodeId() = default;

    friend std::ostream& operator<<(std::ostream& os, NodeId const& o) {
        os << o.row() << ":" << o.col();
        if (o.code_ & NODE_ATTR_MASK) {
            os << "+";
        }
        return os;
    }
};

#endif  // NODE_ID_HPP
