#ifndef NODE_BDD_REDUCER_HPP
#define NODE_BDD_REDUCER_HPP

#include <cassert>             // for assert
#include <cstddef>             // for size_t
#include <ext/alloc_traits.h>  // for __alloc_traits<>::value_type
#include <memory>              // for allocator_traits<>::value_type
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/reverse.hpp>
#include <span>                  // for span
#include <unordered_set>         // for unordered_set
#include <vector>                // for vector
#include "NodeBddTable.hpp"      // for TableHandler, NodeTableEntity
#include "NodeId.hpp"            // for NodeId
#include "util/MyHashTable.hpp"  // for MyHashDefault

template <typename T, bool BDD, bool ZDD>
class DdReducer {
    TableHandler<T>                   oldDiagram;
    NodeTableEntity<T>&               input;
    TableHandler<T>                   newDiagram;
    NodeTableEntity<T>&               output;
    std::vector<std::vector<NodeId>>  newIdTable;
    std::vector<std::vector<NodeId*>> rootPtr;
    size_t                            counter = 1;

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

        for (auto i = 2UL; i < input.numRows(); ++i) {
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

                    for (size_t bb = ZDD ? 1UL : 0UL; bb < 2; ++bb) {
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
    void reduce(size_t i) {
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
    void algorithmR(size_t i) {
        makeReadyForSequentialReduction();
        size_t const m = input[i].size();
        // T* const     tt = input[i].data();

        std::vector<NodeId>& newId = newIdTable[i];
        newId.resize(m);

        // for (size_t j = m - 1; j + 1 > 0; --j) {
        for (auto&& [j, f] :
             input[i] | ranges::views::enumerate | ranges::views::reverse) {
            auto& f0 = f[0];
            auto& f1 = f[1];

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

            auto& g0 = input[i][j][0];
            assert(newId[j].row() == counter + 1);
            newId[j] = NodeId(counter, mm++, g0.hasEmpty());
        }

        auto const& levels = input.lowerLevels(counter);
        for (const auto& t : levels) {
            input[t].clear();
        }

        if (mm > 0U) {
            output.initRow(counter, mm);
            std::span<T> nt{output[counter].data(), output[counter].size()};

            for (size_t j = 0; j < m; ++j) {
                // NodeId const& f0 = tt[j].branch[0];
                auto const& f1 = input[i][j][1];

                if (ZDD && f1 == 0) {  // forwarded
                    assert(newId[j].row() < counter);
                } else {
                    assert(newId[j].row() == counter);
                    auto k = newId[j].col();
                    nt[k] = input[i][j];
                    nt[k].set_node_id_label(newId[j]);
                    if (nt[k].get_ptr_node_id() != 0) {
                        nt[k].set_ptr_node_id(newId[j]);
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
    void algorithmZdd(size_t i) {
        makeReadyForSequentialReduction();
        size_t const       m = input[i].size();
        std::span<T> const tt{input[i].data(), input[i].size()};
        NodeId const       mark(i, m);

        auto& newId = newIdTable[i];
        newId.resize(m);

        for (size_t j = m - 1; j + 1 > 0; --j) {
            auto& f0 = tt[j][0];
            auto& f1 = tt[j][1];

            if (f0.row() != 0) {
                f0 = newIdTable[f0.row()][f0.col()];
            }
            if (f1.row() != 0) {
                f1 = newIdTable[f1.row()][f1.col()];
            }

            if (ZDD && f1 == 0) {
                newId[j] = f0;
            } else {
                auto& f00 = input.child(f0, 0UL);
                auto& f01 = input.child(f0, 1UL);

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
            auto const& levels = input.lowerLevels(i);
            for (auto& t : levels) {
                newIdTable[t].clear();
            }
        }
        size_t mm = 0;

        for (auto j = 0UL; j < m; ++j) {
            NodeId const f(i, j);
            assert(newId[j].row() <= i + 1);
            if (newId[j].row() <= i) {
                continue;
            }

            for (auto k = j; k < m;) {  // for each g in f0-equivalent list
                assert(j <= k);
                NodeId const g(i, k);
                auto&        g0 = tt[k][0];
                auto&        g1 = tt[k][1];
                auto&        g10 = input.child(g1, 0);
                auto&        g11 = input.child(g1, 1);
                assert(g1 != mark);
                assert(newId[k].row() == i + 1);
                auto next = newId[k].col();

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

        auto const& levels = input.lowerLevels(i);
        for (auto& t : levels) {
            input[t].clear();
        }

        if (mm > 0U) {
            output.initRow(i, mm);
            std::span<T> nt{output[i].data(), output[i].size()};

            for (size_t j = 0; j < m; ++j) {
                auto const& f0 = tt[j][0];
                auto const& f1 = tt[j][1];

                if (f1 == mark) {  // forwarded
                    assert(f0.row() == i);
                    assert(newId[j] == 0);
                    newId[j] = newId[f0.col()];
                } else if (ZDD && f1 == 0) {  // forwarded
                    assert(newId[j].row() < i);
                } else {
                    assert(newId[j].row() == i);
                    auto k = newId[j].col();
                    nt[k] = tt[j];
                    nt[k].set_node_id_label(newId[j]);
                }
            }

            counter++;
        }

        for (auto& k : rootPtr[i]) {
            auto& root = *k;
            root = newId[root.col()];
        }
    }

    /**
     * Reduces one level.
     * @param i level.
     */
    void reduce_(int i) {
        auto const m = input[i].size();
        newIdTable[i].resize(m);
        auto jj = 0UL;

        {
            std::unordered_set<T const*> uniq(m * 2, MyHashDefault<T const*>(),
                                              MyHashDefault<T const*>());

            for (auto j = 0UL; j < m; ++j) {
                auto* const p0 = input[i].data();
                auto&       f = input[i][j];

                // make f canonical
                auto& f0 = f[0];
                f0 = newIdTable[f0.row()][f0.col()];
                auto deletable = BDD ? f0 : 0;
                auto del = BDD || ZDD || (f0 == 0);
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
                    auto const* pp = uniq.add(&f);

                    if (pp == &f) {
                        newIdTable[i][j] = NodeId(i, jj++, f0.hasEmpty());
                    } else {
                        newIdTable[i][j] = newIdTable[i][pp - p0];
                    }
                }
            }
        }

        std::vector<int> const& levels = input.lowerLevels(i);
        for (const auto& t : levels) {
            newIdTable[t].clear();
        }

        output.initRow(i, jj);

        for (auto j = 0UL; j < m; ++j) {
            auto const& ff = newIdTable[i][j];
            if (ff.row() == i) {
                output[i][ff.col()] = input[i][j];
                output[i][ff.col()].child[0] =
                    &(output.node(output[i][ff.col()][0]));
                output[i][ff.col()].child[1] =
                    &(output.node(output[i][ff.col()][1]));
            }
        }

        input[i].clear();

        for (auto k = 0UL; k < rootPtr[i].size(); ++k) {
            NodeId& root = *rootPtr[i][k];
            root = newIdTable[i][root.col()];
        }
    }
};

#endif  // NODE_BDD_REDUCER_HPP
