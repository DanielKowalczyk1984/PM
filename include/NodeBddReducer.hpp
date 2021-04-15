#ifndef NODE_BDD_REDUCER_HPP
#define NODE_BDD_REDUCER_HPP

#include <fmt/core.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <span>
#include <stdexcept>

#include <NodeBdd.hpp>
#include <NodeBddTable.hpp>
#include <unordered_set>
#include <vector>
#include "util/MyHashTable.hpp"
// #include "tdzdd/util/MyList.hpp"
#include "util/MyVector.hpp"

template <typename T, bool BDD, bool ZDD>
class DdReducer {
    TableHandler<T>                   oldDiagram;
    NodeTableEntity<T>&               input;
    TableHandler<T>                   newDiagram;
    NodeTableEntity<T>&               output;
    std::vector<std::vector<NodeId>>  newIdTable;
    std::vector<std::vector<NodeId*>> rootPtr;
    int                               counter = 1;

    struct ReduceNodeInfo {
        NodeBdd<T> children;
        size_t     column;

        [[nodiscard]] size_t hash() const { return children.hash(); }

        bool operator==(ReduceNodeInfo const& o) const {
            return children == o.children;
        }

        friend std::ostream& operator<<(std::ostream&         os,
                                        ReduceNodeInfo const& o) {
            return os << "(" << o.children << " -> " << o.column << ")";
        }
    };

    bool readyForSequentialReduction;

   public:
    explicit DdReducer(TableHandler<T>& diagram, bool useMP = false)
        : oldDiagram(std::move(diagram)),
          input(*oldDiagram),
          newDiagram(input.numRows()),
          output(*newDiagram),
          newIdTable(input.numRows()),
          rootPtr(input.numRows()),
          readyForSequentialReduction(false) {
        diagram = std::move(newDiagram);

        input.initTerminals();
        input.makeIndex(useMP);

        newIdTable[0].resize(2);
        newIdTable[0][0] = 0;
        newIdTable[0][1] = 1;
    }

   private:
    /**
     * Applies the node deletion rules.
     * It is required before serial reduction (Algorithm-R)
     * in order to make lower-level index safe.
     */
    void makeReadyForSequentialReduction() {
        if (readyForSequentialReduction) {
            return;
        }

        for (int i = 2; i < input.numRows(); ++i) {
            // size_t const    m = input[i].size();
            // std::span const tt{input[i].data(), input[i].size()};
            // T* const     tt = input[i].data();

            for (auto& it : input[i]) {
                for (auto& f : it) {
                    // NodeId& f = it[b];
                    if (f.row() == 0) {
                        continue;
                    }

                    NodeId f0 = input.child(f, 0);
                    NodeId deletable = 0;
                    bool   del = true;

                    for (int bb = ZDD ? 1 : 0; bb < 2; ++bb) {
                        if (input.child(f, bb) != deletable) {
                            del = false;
                        }
                    }

                    if (del) {
                        f = f0;
                    }
                }
            }
        }

        input.makeIndex();
        readyForSequentialReduction = true;
    }

   public:
    /**
     * Sets a root node.
     * @param root reference to a root node ID storage.
     */
    void setRoot(NodeId& root) { rootPtr[root.row()].push_back(&root); }

    /**
     * Reduces one level.
     * @param i level.
     * @param useMP use an algorithm for multiple processors.
     */
    void reduce(int i) {
        if (BDD) {
            algorithmZdd(i);
        } else if (ZDD) {
            algorithmR(i);
        }
    }

   private:
    /**
     * Reduces one level using Algorithm-R.
     * @param i level.
     */
    void algorithmR(int i) {
        makeReadyForSequentialReduction();
        size_t const m = input[i].size();
        // T* const     tt = input[i].data();

        std::vector<NodeId>& newId = newIdTable[i];
        newId.resize(m);

        for (size_t j = m - 1; j + 1 > 0; --j) {
            NodeId& f0 = input[i][j][0];
            NodeId& f1 = input[i][j][1];

            if (f0.row() != 0) {
                f0 = newIdTable[f0.row()][f0.col()];
            }
            if (f1.row() != 0) {
                f1 = newIdTable[f1.row()][f1.col()];
            }

            if (ZDD && f1 == 0) {
                newId[j] = f0;
            } else {
                newId[j] =
                    NodeId(counter + 1, m);  // tail of f0-equivalent list
            }
        }

        {
            auto const& levels = input.lowerLevels(counter);
            for (auto& t : levels) {
                newIdTable[t].clear();
            }
        }
        size_t mm = 0;

        for (size_t j = 0; j < m; ++j) {
            assert(newId[j].row() <= counter + 1);
            if (newId[j].row() <= counter) {
                continue;
            }

            NodeId& g0 = input[i][j][0];
            assert(newId[j].row() == counter + 1);
            newId[j] = NodeId(counter, mm++, g0.hasEmpty());
        }

        std::vector<int> const& levels = input.lowerLevels(counter);
        for (auto& t : levels) {
            input[t].clear();
        }

        if (mm > 0u) {
            output.initRow(counter, mm);
            std::span<T> nt{output[counter].data(), output[counter].size()};

            for (size_t j = 0; j < m; ++j) {
                // NodeId const& f0 = tt[j].branch[0];
                NodeId const& f1 = input[i][j][1];

                if (ZDD && f1 == 0) {  // forwarded
                    assert(newId[j].row() < counter);
                } else {
                    assert(newId[j].row() == counter);
                    size_t k = newId[j].col();
                    nt[k] = input[i][j];
                    nt[k].set_node_id_label(newId[j]);
                    if (nt[k].ptr_node_id != nullptr) {
                        *(nt[k].ptr_node_id) = newId[j];
                    }
                }
            }

            counter++;
        }

        for (auto& k : rootPtr[i]) {
            NodeId& root = *k;
            root = newId[root.col()];
        }
    }

    /**
     * Reduces one level using Algorithm-R.
     * @param i level.
     */
    void algorithmZdd(int i) {
        makeReadyForSequentialReduction();
        size_t const       m = input[i].size();
        std::span<T> const tt{input[i].data(), input[i].size()};
        NodeId const       mark(i, m);

        auto& newId = newIdTable[i];
        newId.resize(m);

        for (size_t j = m - 1; j + 1 > 0; --j) {
            NodeId& f0 = tt[j][0];
            NodeId& f1 = tt[j][1];

            if (f0.row() != 0) {
                f0 = newIdTable[f0.row()][f0.col()];
            }
            if (f1.row() != 0) {
                f1 = newIdTable[f1.row()][f1.col()];
            }

            if (ZDD && f1 == 0) {
                newId[j] = f0;
            } else {
                NodeId& f00 = input.child(f0, 0);
                NodeId& f01 = input.child(f0, 1);

                if (f01 != mark) {  // the first touch from this level
                    f01 = mark;     // mark f0 as touched
                    newId[j] = NodeId(i + 1, m);  // tail of f0-equivalent list
                } else {
                    newId[j] = f00;  // next of f0-equivalent list
                }
                f00 = NodeId(i + 1, j);  // new head of f0-equivalent list
            }
        }

        {
            std::vector<int> const& levels = input.lowerLevels(i);
            for (auto& t : levels) {
                newIdTable[t].clear();
            }
        }
        size_t mm = 0;

        for (size_t j = 0; j < m; ++j) {
            NodeId const f(i, j);
            assert(newId[j].row() <= i + 1);
            if (newId[j].row() <= i) {
                continue;
            }

            for (size_t k = j; k < m;) {  // for each g in f0-equivalent list
                assert(j <= k);
                NodeId const g(i, k);
                NodeId&      g0 = tt[k][0];
                NodeId&      g1 = tt[k][1];
                NodeId&      g10 = input.child(g1, 0);
                NodeId&      g11 = input.child(g1, 1);
                assert(g1 != mark);
                assert(newId[k].row() == i + 1);
                size_t next = newId[k].col();

                if (g11 != f) {  // the first touch to g1 in f0-equivalent list
                    g11 = f;     // mark g1 as touched
                    g10 = g;     // record g as a canonical node for <f0,g1>
                    newId[k] = NodeId(i, mm++, g0.hasEmpty());
                } else {
                    g0 = g10;   // make a forward link
                    g1 = mark;  // mark g as forwarded
                    newId[k] = 0;
                }

                k = next;
            }
        }

        std::vector<int> const& levels = input.lowerLevels(i);
        for (auto& t : levels) {
            input[t].clear();
        }

        if (mm > 0u) {
            output.initRow(i, mm);
            std::span<T> nt{output[i].data(), output[i].size()};

            for (size_t j = 0; j < m; ++j) {
                NodeId const& f0 = tt[j][0];
                NodeId const& f1 = tt[j][1];

                if (f1 == mark) {  // forwarded
                    assert(f0.row() == i);
                    assert(newId[j] == 0);
                    newId[j] = newId[f0.col()];
                } else if (ZDD && f1 == 0) {  // forwarded
                    assert(newId[j].row() < i);
                } else {
                    assert(newId[j].row() == i);
                    size_t k = newId[j].col();
                    nt[k] = tt[j];
                    fmt::print("test {}\n", nt[k].coeff_list[1].size());
                    nt[k].set_node_id_label(newId[j]);
                }
            }

            counter++;
        }

        for (auto& k : rootPtr[i]) {
            NodeId& root = *k;
            root = newId[root.col()];
        }
    }

    /**
     * Reduces one level.
     * @param i level.
     */
    void reduce_(int i) {
        size_t const m = input[i].size();
        newIdTable[i].resize(m);
        size_t jj = 0;

        {
            // MyList<ReducNodeInfo> rni;
            // MyHashTable<ReducNodeInfo const*> uniq(m * 2);
            std::unordered_set<NodeBdd<T> const*> uniq(
                m * 2, MyHashDefault<NodeBdd<T> const*>(),
                MyHashDefault<NodeBdd<T> const*>());

            for (size_t j = 0; j < m; ++j) {
                NodeBdd<T>* const p0 = input[i].data();
                NodeBdd<T>&       f = input[i][j];

                // make f canonical
                NodeId& f0 = f[0];
                f0 = newIdTable[f0.row()][f0.col()];
                NodeId deletable = BDD ? f0 : 0;
                bool   del = BDD || ZDD || (f0 == 0);
                for (int b = 1; b < 2; ++b) {
                    NodeId& ff = f[b];
                    ff = newIdTable[ff.row()][ff.col()];
                    if (ff != deletable) {
                        del = false;
                    }
                }

                if (del) {  // f is redundant
                    newIdTable[i][j] = f0;
                } else {
                    NodeBdd<T> const* pp = uniq.add(&f);

                    if (pp == &f) {
                        newIdTable[i][j] = NodeId(i, jj++, f0.hasEmpty());
                    } else {
                        newIdTable[i][j] = newIdTable[i][pp - p0];
                    }
                }
            }
        }

        std::vector<int> const& levels = input.lowerLevels(i);
        for (auto& t : levels) {
            newIdTable[t].clear();
        }

        output.initRow(i, jj);

        for (size_t j = 0; j < m; ++j) {
            NodeId const& ff = newIdTable[i][j];
            if (ff.row() == i) {
                output[i][ff.col()] = input[i][j];
                output[i][ff.col()].child[0] =
                    &(output.node(output[i][ff.col()][0]));
                output[i][ff.col()].child[1] =
                    &(output.node(output[i][ff.col()][1]));
            }
        }

        input[i].clear();

        for (size_t k = 0; k < rootPtr[i].size(); ++k) {
            NodeId& root = *rootPtr[i][k];
            root = newIdTable[i][root.col()];
        }
    }
};

#endif  // NODE_BDD_REDUCER_HPP
