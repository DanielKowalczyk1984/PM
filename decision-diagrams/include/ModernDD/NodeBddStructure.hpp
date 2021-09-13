#ifndef __NODEBDDSTRUCTURE_H__
#define __NODEBDDSTRUCTURE_H__

/*
 * TdZdd: a Top-down/Breadth-first Decision Diagram Manipulation Framework
 * by Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
 * Copyright (c) 2014 ERATO MINATO Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include <array>                                 // for array, array<>::valu...
#include <cassert>                               // for assert
#include <cstddef>                               // for size_t
#include <ext/alloc_traits.h>                    // for __alloc_traits<>::va...
#include <memory>                                // for allocator_traits<>::...
#include <range/v3/iterator/basic_iterator.hpp>  // for basic_iterator, oper...
#include <range/v3/view/drop.hpp>                // for drop, drop_fn
#include <range/v3/view/filter.hpp>              // for filter
#include <range/v3/view/iota.hpp>                // for iota_view<>::cursor
#include <range/v3/view/join.hpp>                // for join
#include <range/v3/view/reverse.hpp>             // for reverse
#include <range/v3/view/take.hpp>                // for take, take_fn
#include <set>                                   // for set
#include <utility>                               // for move, pair
#include <vector>                                // for vector
#include "NodeBase.hpp"                          // for InitializedNode
#include "NodeBddBuilder.hpp"                    // for DdBuilder, ZddSubsetter
#include "NodeBddEval.hpp"                       // for Eval
#include "NodeBddReducer.hpp"                    // for DdReducer
#include "NodeBddSpec.hpp"                       // for DdSpec, DdSpecBase
#include "NodeBddTable.hpp"                      // for TableHandler
#include "NodeId.hpp"                            // for NodeId
#include "util/DataTable.hpp"                    // for DataTable
#include "util/MyHashTable.hpp"                  // for MyHashMap

/**
 * Ordered n-ary decision diagram structure.
 * @tparam ARITY arity of the nodes.
 */
template <typename T>
class DdStructure : public DdSpec<DdStructure<T>, NodeId> {
    TableHandler<T> diagram;  ///< The diagram structure.
    NodeId          root_{};  ///< Root node ID.

   public:
    /**
     * Default constructor.
     */
    DdStructure() = default;
    DdStructure(const DdStructure<T>&) = default;
    DdStructure(DdStructure<T>&&) noexcept = default;
    DdStructure<T>& operator=(const DdStructure<T>&) = delete;
    DdStructure<T>& operator=(DdStructure<T>&&) noexcept = default;
    ~DdStructure() = default;

    /**
     * Universal ZDD constructor.
     * @param n the number of variables.
     */
    explicit DdStructure(int n) : diagram(n + 1), root_(1) {
        assert(n >= 0);
        auto&  table = *diagram;
        NodeId f(1);

        for (auto i = 1UL; i <= n; ++i) {
            table.initRow(i, 1);
            for (auto& it : table[i][0]) {
                it = f;
            }
            f = NodeId(i, 0);
        }

        root_ = f;
    }

    /**
     * DD construction.
     * @param spec DD spec.
     */
    template <typename SPEC>
    explicit DdStructure(DdSpecBase<SPEC> const& spec) {
        construct_(spec.entity());
    }

   private:
    template <typename SPEC>
    void construct_(SPEC const& spec) {
        DdBuilder<SPEC, T> zc(spec, diagram);
        int                n = zc.initialize(root_);

        if (n > 0) {
            for (auto i = size_t(n); i > 0UL; --i) {
                zc.construct(i);
            }
        }
    }

   public:
    /**
     * ZDD subsetting.
     * @param spec ZDD spec.
     */
    template <typename SPEC>
    void zddSubset(SPEC const& spec) {
#ifdef _OPENMP
        if (useMP)
            zddSubsetMP_(spec.entity());
        else
#endif
            zddSubset_(spec.entity());
    }

   private:
    template <typename SPEC>
    void zddSubset_(SPEC const& spec) {
        TableHandler<T>       tmpTable;
        ZddSubsetter<T, SPEC> zs(diagram, spec, tmpTable);
        int                   n = zs.initialize(root_);

        if (n > 0) {
            for (int i = n; i > 0; --i) {
                zs.subset(i);
                // diagram.derefLevel(i);
            }
        }

        diagram = std::move(tmpTable);
    }

   public:
    /**
     * Gets the root node.
     * @return root node ID.
     */
    NodeId& root() { return root_; }

    /**
     * Gets the root node.
     * @return root node ID.
     */
    [[nodiscard]] NodeId root() const { return root_; }

    /**
     * Gets a child node.
     * @param f parent node ID.
     * @param b branch number.
     * @return child node ID.
     */
    [[nodiscard]] NodeId child(NodeId f, size_t b) const {
        return diagram->child(f, b);
    }

    /**
     * Gets the diagram.
     * @return the node table handler.
     */
    TableHandler<T>& getDiagram() { return diagram; }

    /**
     * Gets the diagram.
     * @return the node table handler.
     */
    TableHandler<T> const& getDiagram() const { return diagram; }

    /**
     * Gets the level of the root node.
     * @return the level of root ZDD variable.
     */
    [[nodiscard]] size_t topLevel() const { return root_.row(); }

    /**
     * Gets the number of nonterminal nodes.
     * @return the number of nonterminal nodes.
     */
    [[nodiscard]] size_t size() const { return diagram->size(); }

    /**
     * Checks if DD is a 0-terminal only.
     * @return true if DD is a 0-terminal only.
     */
    [[nodiscard]] bool empty() const { return root_ == 0; }

    /**
     * Checks structural equivalence with another DD.
     * @return true if they have the same structure.
     */
    bool operator==(DdStructure const& o) const {
        auto n = root_.row();
        if (n != o.root_.row()) {
            return false;
        }
        if (n == 0) {
            return root_ == o.root_;
        }
        if (root_ == o.root_ && &*diagram == &*o.diagram) {
            return true;
        }
        if (size() > o.size()) {
            return o.operator==(*this);
        }

        MyHashMap<NodeBase, size_t> uniq;
        DataTable<NodeId>           equiv(n + 1);
        {
            size_t om = (*o.diagram)[0].size();
            equiv[0].resize(om);
            for (size_t j = 0; j < om; ++j) {
                equiv[0][j] = j;
            }
        }

        for (auto i = 1UL; i <= n; ++i) {
            auto m = (*diagram)[i].size();
            uniq.initialize(m * 2);

            for (auto j = 0; j < m; ++j) {
                uniq[(*diagram)[i][j]] = j;
            }

            auto om = (*o.diagram)[i].size();
            equiv[i].resize(om);

            for (auto j = 0UL; j < om; ++j) {
                NodeBase node;

                for (auto b = 0UL; b < 2; ++b) {
                    NodeId f = (*o.diagram)[i][j][b];
                    node[b] = equiv[f.row()][f.col()];
                }

                auto* p = uniq.getValue(node);
                equiv[i][j] = NodeId(i, (p != nullptr) ? *p : m);
            }
        }

        return root_ == equiv[o.root_.row()][o.root_.col()];
    }

    /**
     * Checks structural in-equivalence with another DD.
     * @return true if they have the different structure.
     */
    bool operator!=(DdStructure const& o) const { return !operator==(o); }

    /**
     * QDD reduction.
     * No node deletion rule is applied.
     */
    void qddReduce() { reduce<false, false>(); }

    /**
     * BDD reduction.
     * The node of which two edges points to the identical node is deleted.
     */
    void bddReduce() { reduce<true, false>(); }

    /**
     * ZDD reduction.
     * The node of which 1-edge points to the 0-terminal is deleted.
     */
    void compressBdd() { reduce<false, true>(); }

    void reduceZdd() { reduce<true, true>(); }

    /**
     * BDD/ZDD reduction.
     * @tparam BDD enable BDD reduction.
     * @tparam ZDD enable ZDD reduction.
     */
    template <bool BDD, bool ZDD>
    void reduce() {
        auto n = root_.row();

        DdReducer<T, BDD, ZDD> zr(diagram);
        zr.setRoot(root_);

        for (auto i : ranges::views::ints(1UL, n + 1)) {
            zr.reduce(i);
        }
    }

    // /**
    //  * Transforms a BDD into a ZDD.
    //  * @param numVars the number of variables.
    //  */
    // DdStructure bdd2zdd(int numVars) const {
    //     return DdStructure(
    //             ZddLookahead<BddUnreduction<DdStructure> >(
    //                     BddUnreduction<DdStructure>(*this, numVars)), useMP);
    // }

    // /**
    //  * Transforms a ZDD into a BDD.
    //  * @param numVars the number of variables.
    //  */
    // DdStructure zdd2bdd(int numVars) const {
    //     return DdStructure(
    //             BddLookahead<ZddUnreduction<DdStructure> >(
    //                     ZddUnreduction<DdStructure>(*this, numVars)), useMP);
    // }

    // /**
    //  * Counts the number of minterms of the function represented by this BDD.
    //  * @param numVars the number of input variables of the function.
    //  * @return the number of itemsets.
    //  */
    // std::string bddCardinality(int numVars) const {
    //     return evaluate(BddCardinality<std::string,ARITY>(numVars));
    // }

    // /**
    //  * Counts the number of sets in the family of sets represented by this
    //  ZDD.
    //  * @return the number of itemsets.
    //  */
    // std::string zddCardinality() const {
    //     return evaluate(ZddCardinality<std::string,ARITY>());
    // }

    template <typename R>
    R evaluate_backward(Eval<T, R>& evaluator) {
        auto  n = root_.row();
        auto& work = *diagram;
        evaluator.set_table(&(*diagram));

        if (this->size() == 0) {
            // fmt::print("empty DDstructure\n");
            R retval{};
            return retval;
        }

        evaluator.initializerootnode(work.node(1));
        for (auto& it : work | ranges::views::take(n + 1) |
                            ranges::views::drop(1) | ranges::views::join) {
            evaluator.initializenode(it);
            evaluator.evalNode(it);
        }

        return evaluator.get_objective(work.node(root()));
    }

    template <typename R>
    void compute_labels_backward(Eval<T, R>& evaluator) {
        auto  n = root_.row();
        auto& work = *(getDiagram());
        evaluator.set_table(&(*diagram));

        if (this->size() == 0) {
            // fmt::print("empty DDstructure\n");
            return;
        }

        evaluator.initializerootnode(work.node(1));
        // for (int i = 1; i <= n; ++i) {
        for (auto& it : work | ranges::views::take(n + 1) |
                            ranges::views::drop(1) | ranges::views::join) {
            evaluator.initializenode(it);
            evaluator.evalNode(it);
        }
        // }
    }

    template <typename R>
    R evaluate_forward(Eval<T, R>& evaluator) {
        auto  n = root_.row();
        auto& work = *(getDiagram());
        evaluator.set_table(&(*diagram));

        if (this->size() == 0) {
            // fmt::print("empty DDstructure\n");
            R retval;
            return retval;
        }

        /**
         * Initialize nodes of the DD
         */
        evaluator.initializerootnode(work.node(root()));
        // for (int i = n - 1; i >= 0; i--) {
        for (auto& it : work | ranges::views::take(n) | ranges::views::reverse |
                            ranges::views::join) {
            evaluator.initializenode(it);
        }
        // }

        /**
         * Compute all the node of DD
         */
        // for (int i = n; i > 0; i--) {
        for (auto& it : work | ranges::views::take(n + 1) |
                            ranges::views::drop(1) | ranges::views::reverse |
                            ranges::views::join) {
            evaluator.evalNode(it);
        }
        // }

        /**
         * Return the optimal solution
         */
        return evaluator.get_objective(work.node(1));
    }

    template <typename R>
    void compute_labels_forward(Eval<T, R>& evaluator) {
        auto  n = root_.row();
        auto& work = *(getDiagram());
        evaluator.set_table(&(*diagram));

        if (this->size() == 0) {
            // fmt::print("empty DDstructure\n");
            return;
        }

        /**
         * Initialize nodes of the DD
         */
        evaluator.initializerootnode(work.node(root()));
        for (auto& it : work | ranges::views::take(n) | ranges::views::reverse |
                            ranges::views::join) {
            evaluator.initializenode(it);
        }

        /**
         * Compute all the node of DD
         */
        // for (int i = n; i > 0; i--) {
        for (auto& it : work | ranges::views::take(n + 1) |
                            ranges::views::drop(1) | ranges::views::reverse |
                            ranges::views::join) {
            evaluator.evalNode(it);
        }
        // }
    }

    /**
     * Iterator on a set of integer vectors represented by a DD.
     */
    class const_iterator {
        struct Selection {
           private:
            NodeId node{};
            bool   val{false};

           public:
            Selection() = default;

            Selection(NodeId _node, bool _val) : node(_node), val(_val) {}

            bool operator==(Selection const& o) const {
                return node == o.node && val == o.val;
            }
        };

        DdStructure const&                   dd;
        int                                  cursor;
        std::vector<std::pair<NodeId, bool>> path;
        std::set<int>                        itemset;

       public:
        const_iterator(DdStructure const& _dd, bool begin)
            : dd(_dd),
              cursor(begin ? -1 : -2),
              path(),
              itemset() {
            if (begin) {
                next(dd.root_);
            }
        }

        const_iterator& operator++() {
            next(NodeId(0, 0));
            return *this;
        }

        std::set<int> const& operator*() const { return itemset; }

        std::set<int> const* operator->() const { return &itemset; }

        bool operator==(const_iterator const& o) const {
            return cursor == o.cursor && path == o.path;
        }

        bool operator!=(const_iterator const& o) const {
            return !operator==(o);
        }

       private:
        void next(NodeId f) {
            itemset.clear();

            for (;;) {
                while (f > 1) { /* down */
                    auto const& s = (*dd.diagram)[f.row()][f.col()];

                    if (s[0] != 0) {
                        cursor = path.size();
                        path.emplace_back(f, false);
                        f = s[0];
                    } else {
                        path.emplace_back(f, true);
                        f = s[1];
                    }
                }

                if (f == 1) {
                    break; /* found */
                }

                for (; cursor >= 0; --cursor) { /* up */
                    auto&       sel = path[cursor];
                    auto const& ss =
                        (*dd.diagram)[sel.first.row()][sel.first.col()];
                    if (sel.second == false && ss[1] != 0) {
                        f = sel.first;
                        sel.second = true;
                        path.resize(cursor + 1);
                        f = dd.diagram->child(f, 1);
                        break;
                    }
                }

                if (cursor < 0) { /* end() state */
                    cursor = -2;
                    path.clear();
                    return;
                }
            }

            for (auto& it : path | ranges::views::filter([](const auto& tmp) {
                                return tmp.second;
                            })) {
                itemset.insert(it.first.row());
            }
        }
    };

    /**
     * Returns an iterator to the first instance,
     * which is viewed as a collection of item numbers.
     * Supports binary ZDDs only.
     * @return iterator to the first instance.
     */
    const_iterator begin() const { return const_iterator(*this, true); }

    /**
     * Returns an iterator to the element following the last instance.
     * Supports binary ZDDs only.
     * @return iterator to the instance following the last instance.
     */
    const_iterator end() const { return const_iterator(*this, false); }

    /**
     * Implements DdSpec.
     */
    size_t getRoot(NodeId& f) const {
        f = root_;
        return (f == 1) ? size_t(-1) : f.row();
    }

    /**
     * Implements DdSpec.
     */
    size_t getChild(NodeId& f, [[maybe_unused]] int level, size_t value) const {
        assert(level > 0 && level == f.row());
        assert(0 <= value && value < 2);
        f = child(f, value);
        return (f.row() > 0) ? f.row() : -f.col();
    }

    /**
     * Implements DdSpec.
     */
    [[nodiscard]] size_t hashCode(NodeId const& f) const { return f.hash(); }
};

#endif  // __NODEBDDSTRUCTURE_H__