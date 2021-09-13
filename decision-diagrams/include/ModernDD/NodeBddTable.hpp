#ifndef NODE_BDD_TABLE_HPP
#define NODE_BDD_TABLE_HPP

#include <cassert>             // for assert
#include <cstddef>             // for size_t
#include <memory>              // for allocator, allocator_traits<>::value_type
#include <ostream>             // for operator<<, ostream, basic_ostream
#include <stdexcept>           // for runtime_error
#include <string>              // for operator<<, char_traits, string
#include <vector>              // for vector, _Bit_reference, vector<>::refe...
#include "NodeId.hpp"          // for NodeId, operator<<
#include "util/DataTable.hpp"  // for DataTable

template <typename T>
using data_table_node = DataTable<T>;
template <typename T>
using my_vector = std::vector<T>;

template <typename T>
class NodeTableEntity : public data_table_node<T> {
    mutable my_vector<my_vector<size_t>> higherLevelTable;
    mutable my_vector<my_vector<size_t>> lowerLevelTable;

   public:
    /**
     * Constructor.
     * @param n the number of rows.
     */
    // explicit NodeTableEntity(int n = 1) : data_table_node<T>(n) {
    //     assert(n >= 1);
    //     initTerminals();
    // }

    explicit NodeTableEntity(size_t n = 1) : data_table_node<T>(n) {
        assert(n >= 1);
        initTerminals();
    }
    NodeTableEntity(const NodeTableEntity<T>&) = default;
    NodeTableEntity<T>& operator=(const NodeTableEntity<T>&) = delete;
    NodeTableEntity<T>& operator=(NodeTableEntity<T>&&) noexcept = default;
    NodeTableEntity(NodeTableEntity<T>&&) noexcept = default;
    ~NodeTableEntity() = default;

    /**
     * Clears and initializes the table.
     * @param n the number of rows.
     */
    void init(int n) {
        assert(n >= 1);
        data_table_node<T>::init(n);
        initTerminals();
    }

    /**
     * Initializes the terminal nodes.
     */
    void initTerminals() {
        auto& t = (*this)[0];
        t.resize(2);

        for (auto j = 0UL; j < 2; ++j) {
            t[j] = T(j, j);
        }
    }

    //    /**
    //     * Gets the variable ID at a given level.
    //     * @param level level.
    //     * @return variable ID.
    //     */
    //    int varAtLevel(int level) const {
    //        assert(0 <= level && level <= numVars());
    //        return (level == 0) ? INT_MAX : numVars() - level;
    //    }
    //
    //    /**
    //     * Gets the level of a variable.
    //     * @param var variable ID.
    //     * @return level.
    //     */
    //    int levelOfVar(int var) const {
    //        assert((0 <= var && var < numVars()) || var == INT_MAX);
    //        return (var == INT_MAX) ? 0 : numVars() - var;
    //    }

    /**
     * Gets the number of nonterminal nodes.
     * @return the number of nonterminal nodes.
     */
    size_t size() const { return this->totalSize() - (*this)[0].size(); }

    /**
     * Gets the number of ZDD variables.
     * @return the number of ZDD variables.
     */
    int numVars() const { return this->numRows() - 1; }

    /**
     * Changes the number of ZDD variables
     * by shifting up/down the levels of existing variables.
     * @param n required number of variables.
     */
    void stretchBottom(int n) {
        int n0 = numVars();
        int d = n - n0;

        if (d > 0) {
            this->setNumRows(n + 1);

            for (int i = n0; i > 0; --i) {
                size_t m = (*this)[i].size();
                this->initRow(i + d, m);

                for (size_t j = 0; j < m; ++j) {
                    for (int b = 0; b < 2; ++b) {
                        auto ff = child(i, j, b);
                        auto ii = ff.row();
                        child(i + d, j, b) =
                            (ii == 0UL) ? ff : NodeId(ii + d, ff.col());
                    }
                }

                this->initRow(i, 0);
            }
        } else if (d < 0) {
            for (int i = 1 - d; i <= n0; ++i) {
                size_t m = (*this)[i].size();
                this->initRow(i + d, m);

                for (size_t j = 0; j < m; ++j) {
                    for (int b = 0; b < 2; ++b) {
                        NodeId ff = child(i, j, b);
                        int    ii = ff.row();
                        child(i + d, j, b) = (ii == 0) ? ff
                                             : (ii + d <= 0)
                                                 ? 1
                                                 : NodeId(ii + d, ff.col());
                    }
                }

                this->initRow(i, 0);
            }

            this->setNumRows(n + 1);
        }
    }

    /**
     * Gets a node.
     * @param f node ID.
     * @return node @p f.
     */
    T const& node(NodeId f) const { return (*this)[f.row()][f.col()]; }

    /**
     * Gets a reference to a node.
     * @param f node ID.
     * @return node @p f.
     */
    T& node(NodeId f) { return (*this)[f.row()][f.col()]; }

    T* node_ptr(NodeId f) { return &(*this)[f.row()][f.col()]; }

    /**
     * Gets a child node ID.
     * @param f parent node ID.
     * @param b child branch.
     * @return the @p b-child of @p f.
     */
    NodeId child(NodeId f, size_t b) const {
        return child(f.row(), f.col(), b);
    }

    /**
     * Gets a reference to a child node ID.
     * @param f parent node ID.
     * @param b child branch.
     * @return the @p b-child of @p f.
     */
    NodeId& child(NodeId f, size_t b) { return child(f.row(), f.col(), b); }

    /**
     * Gets a child node ID.
     * @param i parent row.
     * @param j parent column.
     * @param b child branch.
     * @return the @p b-child of the parent.
     */
    NodeId child(size_t i, size_t j, size_t b) const {
        assert(b < 2UL);
        return (*this)[i][j][b];
    }

    /**
     * Gets a reference to a child node ID.
     * @param i parent row.
     * @param j parent column.
     * @param b child branch.
     * @return the @p b-child of the parent.
     */
    NodeId& child(size_t i, size_t j, size_t b) {
        assert(b < 2UL);
        return (*this)[i][j][b];
    }

    /**
     * Gets a descendant node ID by tracing 0-edges.
     * @param f parent node ID.
     * @param stopLevel level to stop going down.
     * @return reached node ID.
     */
    NodeId getZeroDescendant(NodeId f, int stopLevel) const {
        assert(0 <= stopLevel);

        if (stopLevel == 0 && f.hasEmpty()) {
            return 1;
        }

        while (f.row() > stopLevel) {
            f = child(f, 0);
        }

        return f;
    }

    /**
     * Deletes current index information.
     */
    void deleteIndex() {
        higherLevelTable.clear();
        lowerLevelTable.clear();
    }

    /**
     * Makes index information.
     * @param useMP use an algorithm for multiple processors.
     */
    void makeIndex(bool useMP = false) const {
        size_t const n = this->numRows() - 1;
        higherLevelTable.clear();
        higherLevelTable.resize(n + 1);
        lowerLevelTable.clear();
        lowerLevelTable.resize(n + 1);
        my_vector<bool> lowerMark(n + 1);

        for (auto i = n; i >= 1; --i) {
            auto const&     node = (*this)[i];
            auto const      m = node.size();
            auto            lowest = i;
            my_vector<bool> myLower(n + 1);

            if (useMP) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
                for (intmax_t j = 0; j < intmax_t(m); ++j) {
                    for (int b = 0; b < 2; ++b) {
                        int const ii = node[j][b].row();
                        if (ii == 0)
                            continue;
                        if (ii < lowest) {
#pragma omp critical
                            if (ii < lowest)
                                lowest = ii;
                        }
                        if (!lowerMark[ii]) {
                            myLower[ii] = true;
                            lowerMark[ii] = true;
                        }
                    }
                }
#endif
            } else {
                for (size_t j = 0; j < m; ++j) {
                    for (auto b = 0UL; b < 2; ++b) {
                        auto const ii = node[j][b].row();

                        if (ii == 0UL) {
                            continue;
                        }

                        if (ii < lowest) {
                            lowest = ii;
                        }

                        if (!lowerMark[ii]) {
                            myLower[ii] = true;
                            lowerMark[ii] = true;
                        }
                    }
                }
            }

            higherLevelTable[lowest].push_back(i);
            auto& lower = lowerLevelTable[i];

            for (size_t ii = lowest; ii < i; ++ii) {
                if (myLower[ii]) {
                    lower.push_back(ii);
                }
            }
        }
    }

    /**
     * Returns a collection of the higher levels that directly refers
     * the given level and that does not refer any lower levels.
     * @param level the level.
     */
    my_vector<size_t> const& higherLevels(int level) const {
        if (higherLevelTable.empty()) {
            makeIndex();
        }

        return higherLevelTable[level];
    }

    /**
     * Returns a collection of the lower levels that are referred
     * by the given level and that are not referred directly by
     * any higher levels.
     * @param level the level.
     */
    my_vector<size_t> const& lowerLevels(size_t level) const {
        if (lowerLevelTable.empty()) {
            makeIndex();
        }

        return lowerLevelTable[level];
    }

    /**
     * Dumps the node table in Graphviz (dot) format.
     * @param os output stream.
     * @param title title label.
     */
    void dumpDot(std::ostream& os, const std::string& title = "") const {
        os << "digraph \"" << title << "\" {\n";

        for (int i = this->numRows() - 1; i >= 1; --i) {
            os << "  " << i << " [shape=none];\n";
        }

        for (int i = this->numRows() - 2; i >= 1; --i) {
            os << "  " << (i + 1) << " -> " << i << " [style=invis];\n";
        }

        if (!title.empty()) {
            os << "  labelloc=\"t\";\n";
            os << "  label=\"" << title << "\";\n";
        }

        bool terminal1 = false;

        for (auto i = this->numRows() - 1; i > 0; --i) {
            size_t m = (*this)[i].size();

            for (size_t j = 0; j < m; ++j) {
                NodeId f = NodeId(i, j);
                os << "  \"" << f << "\";\n";

                for (int b = 0; b < 2; ++b) {
                    NodeId ff = child(i, j, b);
                    bool   aa = ff.getAttr();

                    if (ff == 0) {
                        continue;
                    }

                    if (ff == 1) {
                        terminal1 = true;
                        os << "  \"" << f << R"(" -> "$")";
                    } else {
                        ff.setAttr(false);
                        os << "  \"" << f << "\" -> \"" << ff << "\"";
                    }

                    os << " [style=";

                    if (b == 0) {
                        os << "dashed";
                    } else {
                        os << "solid";
                        // if (ARITY > 2) {
                        //     os << ",color="
                        //             << ((b == 1) ? "blue" :
                        //                 (b == 2) ? "red" : "green");
                        // }
                    }

                    if (aa) {
                        os << ",arrowtail=dot";
                    }

                    os << "];\n";
                }
            }

            if (terminal1) {
                os << "  \"$\" [shape=square,label=\"⊤\"];\n";
            }

            os << "  {rank=same; " << i;

            for (size_t j = 0; j < m; ++j) {
                os << "; \"" << NodeId(i, j) << "\"";
            }

            os << "}\n";
        }

        os << "}\n";
        os.flush();
    }
};

template <typename T = double>
class TableHandler {
    struct Object {
       private:
        unsigned           refCount;
        NodeTableEntity<T> entity;

        friend class TableHandler<T>;

       public:
        explicit Object(size_t n) : refCount(1), entity(n) {}

        explicit Object(const NodeTableEntity<T>& _entity)
            : refCount(1),
              entity(_entity) {}

        void ref() {
            ++refCount;
            // fmt::print("We are using ref {}\n", refCount);

            if (refCount == 0) {
                throw std::runtime_error("Too many references");
            }
        }

        void deref() {
            // fmt::print("we are using deref {}\n", refCount);
            --refCount;

            if (refCount == 0) {
                delete this;
            }
        }
    };

    Object* pointer;

   public:
    explicit TableHandler(size_t n = 1) : pointer(new Object(n)) {}

    TableHandler(TableHandler<T> const& o)
        : pointer(new Object(o.pointer->entity)){};

    //     : pointer(o.pointer) {
    //     fmt::print("we are using copy constructor\n");
    //     pointer->ref();
    //     fmt::print("we are using copy constructor end\n");
    // }

    TableHandler<T>& operator=(TableHandler<T> const& o) = delete;

    //     {
    //     fmt::print("we are using copy assignment\n");
    //     pointer->deref();
    //     pointer = o.pointer;
    //     pointer->ref();
    //     return *this;
    // }

    TableHandler<T>& operator=(TableHandler<T>&& o) noexcept {
        if (this != &o) {
            delete pointer;
            pointer = o.pointer;
            o.pointer = nullptr;
        }
        return *this;
    }

    TableHandler(TableHandler<T>&& o) noexcept : pointer(o.pointer) {
        o.pointer = nullptr;
    }

    ~TableHandler() {
        if (pointer) {
            pointer->deref();
        }
    }

    NodeTableEntity<T>& operator*() const { return pointer->entity; }

    NodeTableEntity<T> const* operator->() const { return &pointer->entity; }

    /**
     * Make the table unshared.
     * @return writable reference to the private table.
     */
    // NodeTableEntity<T>& privateEntity() {
    //     if (pointer->refCount >= 2) {
    //         pointer->deref();
    //         pointer = new Object(pointer->entity);
    //     }

    //     return pointer->entity;
    // }

    /**
     * Clear and initialize the table.
     * @param n the number of rows.
     * @return writable reference to the private table.
     */
    // NodeTableEntity<T>& init(int n = 1) {
    //     if (pointer->refCount == 1) {
    //         pointer->entity.init(n);
    //     } else {
    //         pointer->deref();
    //         pointer = new Object(n);
    //     }

    //     return pointer->entity;
    // }

    /**
     * Clear a row if it is not shared.
     * @param i row index.
     */
    // void derefLevel(int i) {
    //     if (pointer->refCount == 1) {
    //         pointer->entity[i].clear();
    //     }
    // }
};

#endif  // NODE_BDD_TABLE_HPP
