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

#pragma once

#include <algorithm>
#include <cassert>
#include <climits>
#include <ostream>
#include <set>
#include <stdexcept>
#include <vector>

#include <NodeBddTable.hpp>
#include <NodeBddBuilder.hpp>
#include <NodeBddReducer.hpp>

#include "tdzdd/DdEval.hpp"
#include "tdzdd/DdSpec.hpp"
#include "tdzdd/util/MyHashTable.hpp"
#include "tdzdd/util/MyVector.hpp"

/**
 * Ordered n-ary decision diagram structure.
 * @tparam ARITY arity of the nodes.
 */
template<typename T = double>
class DdStructure: public tdzdd::DdSpec<DdStructure<T>,tdzdd::NodeId,2> {
    NodeTableHandler<T> diagram; ///< The diagram structure.
    nodeid root_;                    ///< Root node ID.

public:
    /**
     * Default constructor.
     */
    DdStructure() : root_(0) { }

    /**
     * Universal ZDD constructor.
     * @param n the number of variables.
     * @param useMP use algorithms for multiple processors.
     */
    explicit DdStructure(int n) :
            diagram(n + 1), root_(1) {
        assert(n >= 0);
        NodeTableEntity<T>& table = diagram.privateEntity();
        nodeid f(1);

        for (int i = 1; i <= n; ++i) {
            table.initRow(i, 1);
            for (int b = 0; b < 2; ++b) {
                table[i][0].branch[b] = f;
            }
            f = nodeid(i, 0);
        }

        root_ = f;
    }

    /**
     * DD construction.
     * @param spec DD spec.
     * @param useMP use algorithms for multiple processors.
     */
    template<typename SPEC>
    explicit DdStructure(tdzdd::DdSpecBase<SPEC,2> const& spec) {
        construct_(spec.entity());
    }

private:
    template<typename SPEC>
    void construct_(SPEC const& spec) {
        DdBuilder<SPEC> zc(spec, diagram);
        int n = zc.initialize(root_);

        if (n > 0) {
            for (int i = n; i > 0; --i) {
                zc.construct(i);
            }
        }
    }

public:
    /**
     * ZDD subsetting.
     * @param spec ZDD spec.
     */
//     template<typename SPEC>
//     void zddSubset(DdSpecBase<SPEC,ARITY> const& spec) {
// #ifdef _OPENMP
//         if (useMP) zddSubsetMP_(spec.entity());
//         else
// #endif
//         zddSubset_(spec.entity());
//     }

private:
    template<typename SPEC>
    void zddSubset_(SPEC const& spec) {
        NodeTableHandler<T> tmpTable;
        tdzdd::ZddSubsetter<SPEC> zs(diagram, spec, tmpTable);
        int n = zs.initialize(root_);

        if (n > 0) {
            for (int i = n; i > 0; --i) {
                diagram.derefLevel(i);
            }
        }

        diagram = tmpTable;
    }

public:
    /**
     * Gets the root node.
     * @return root node ID.
     */
    nodeid& root() {
        return root_;
    }

    /**
     * Gets the root node.
     * @return root node ID.
     */
    nodeid root() const {
        return root_;
    }

    /**
     * Gets a child node.
     * @param f parent node ID.
     * @param b branch number.
     * @return child node ID.
     */
    nodeid child(nodeid f, int b) const {
        return diagram->child(f, b);
    }

    /**
     * Gets the diagram.
     * @return the node table handler.
     */
    NodeTableHandler<T>& getDiagram() {
        return diagram;
    }

    /**
     * Gets the diagram.
     * @return the node table handler.
     */
    NodeTableHandler<T> const& getDiagram() const {
        return diagram;
    }

    /**
     * Gets the level of the root node.
     * @return the level of root ZDD variable.
     */
    int topLevel() const {
        return root_.row();
    }

    /**
     * Gets the number of nonterminal nodes.
     * @return the number of nonterminal nodes.
     */
    size_t size() const {
        return diagram->size();
    }

    /**
     * Checks if DD is a 0-terminal only.
     * @return true if DD is a 0-terminal only.
     */
    bool empty() const {
        return root_ == 0;
    }

    /**
     * Checks structural equivalence with another DD.
     * @return true if they have the same structure.
     */
    bool operator==(DdStructure const& o) const {
        int n = root_.row();
        if (n != o.root_.row()) return false;
        if (n == 0) return root_ == o.root_;
        if (root_ == o.root_ && &*diagram == &*o.diagram) return true;
        if (size() > o.size()) return o.operator==(*this);

        tdzdd::MyHashMap<tdzdd::InitializedNode<2>,size_t> uniq;
        tdzdd::DataTable<nodeid> equiv(n + 1);
        {
            size_t om = (*o.diagram)[0].size();
            equiv[0].resize(om);
            for (size_t j = 0; j < om; ++j) {
                equiv[0][j] = j;
            }
        }

        for (int i = 1; i <= n; ++i) {
            size_t m = (*diagram)[i].size();
            uniq.initialize(m * 2);

            for (size_t j = 0; j < m; ++j) {
                uniq[(*diagram)[i][j]] = j;
            }

            size_t om = (*o.diagram)[i].size();
            equiv[i].resize(om);

            for (size_t j = 0; j < om; ++j) {
                tdzdd::InitializedNode<2> node;

                for (int b = 0; b < 2; ++b) {
                    nodeid f = (*o.diagram)[i][j].branch[b];
                    node.branch[b] = equiv[f.row()][f.col()];
                }

                size_t* p = uniq.getValue(node);
                equiv[i][j] = nodeid(i, (p != 0) ? *p : m);
            }
        }

        return root_ == equiv[o.root_.row()][o.root_.col()];
    }

    /**
     * Checks structural inequivalence with another DD.
     * @return true if they have the different structure.
     */
    bool operator!=(DdStructure const& o) const {
        return !operator==(o);
    }

    /**
     * QDD reduction.
     * No node deletion rule is applied.
     */
    void qddReduce() {
        reduce<false,false>();
    }

    /**
     * BDD reduction.
     * The node of which two edges points to the identical node is deleted.
     */
    void bddReduce() {
        reduce<true,false>();
    }

    /**
     * ZDD reduction.
     * The node of which 1-edge points to the 0-terminal is deleted.
     */
    void zddReduce() {
        reduce<false,true>();
    }

    /**
     * BDD/ZDD reduction.
     * @tparam BDD enable BDD reduction.
     * @tparam ZDD enable ZDD reduction.
     */
    template<bool BDD, bool ZDD>
    void reduce() {
        int n = root_.row();

        DdReducer<T,BDD,ZDD> zr(diagram);
        zr.setRoot(root_);

        for (int i = 1; i <= n; ++i) {
            zr.reduce(i);
        }
    }

public:
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
    //  * Counts the number of sets in the family of sets represented by this ZDD.
    //  * @return the number of itemsets.
    //  */
    // std::string zddCardinality() const {
    //     return evaluate(ZddCardinality<std::string,ARITY>());
    // }


    template <typename S, typename R>
    R evaluate_backward(tdzdd::DdEval<S, Node<T>, R> const &evaluator) {
        int            n = root_.row();
        NodeTableHandler<double> &handler_bdd = getDiagram();
        NodeTableEntity<T>& work = handler_bdd.privateEntity();

        if (this->size() == 0) {
            printf("empty DDstructure\n");
            R retval;
            return retval;
        }

        evaluator.initializerootnode(work[0][1]);
        for (int i = 1; i <= n; ++i) {
            for (auto &it : work[i]) {
                evaluator.initializenode(it);
                evaluator.evalNode(it);
            }
        }

        return evaluator.get_objective(work[n][0]);
    }

    template<typename S, typename R>
    R evaluate_forward(tdzdd::DdEval<S, Node<T>, R> const &evaluator) {
        int n = root_.row();
        NodeTableEntity<T> & work = getDiagram().privateEntity();

        if (this->size() == 0){
            printf("empty DDstructure\n");
            R retval;
            return retval;
        }

        /**
         * Initialize nodes of the DD
         */
        evaluator.initializerootnode(work[n][0]);
        for (int i = n - 1; i >= 0; i--) {
            for(auto &it : work[i]){
                evaluator.initializenode(it);
            }
        }

        /**
         * Compute all the node of ZDD
         */
        for (int i = n ; i > 0; i--) {
            for (auto &it: work[i]) {
                evaluator.evalNode(it);
            }
        }

        /**
         * Return the optimal solution
         */
        return evaluator.get_objective(work[0][1]);
    }

    /**
     * Iterator on a set of integer vectors represented by a DD.
     */
    class const_iterator {
        struct Selection {
            nodeid node;
            bool val;

            Selection() :
                    val(false) {
            }

            Selection(nodeid node, bool val) :
                    node(node), val(val) {
            }

            bool operator==(Selection const& o) const {
                return node == o.node && val == o.val;
            }
        };

        DdStructure const& dd;
        int cursor;
        std::vector<Selection> path;
        std::set<int> itemset;

    public:
        const_iterator(DdStructure const& dd, bool begin) :
                dd(dd), cursor(begin ? -1 : -2), path(), itemset() {
            if (begin) next(dd.root_);
        }

        const_iterator& operator++() {
            next(nodeid(0, 0));
            return *this;
        }

        std::set<int> const& operator*() const {
            return itemset;
        }

        std::set<int> const* operator->() const {
            return &itemset;
        }

        bool operator==(const_iterator const& o) const {
            return cursor == o.cursor && path == o.path;
        }

        bool operator!=(const_iterator const& o) const {
            return !operator==(o);
        }

    private:
        void next(nodeid f) {
            itemset.clear();

            for (;;) {
                while (f > 1) { /* down */
                    Node<T> const& s = (*dd.diagram)[f.row()][f.col()];

                    if (s.branch[0] != 0) {
                        cursor = path.size();
                        path.push_back(Selection(f, false));
                        f = s.branch[0];
                    }
                    else {
                        path.push_back(Selection(f, true));
                        f = s.branch[1];
                    }
                }

                if (f == 1) break; /* found */

                for (; cursor >= 0; --cursor) { /* up */
                    Selection& sel = path[cursor];
                    Node<T> const& ss =
                            (*dd.diagram)[sel.node.row()][sel.node.col()];
                    if (sel.val == false && ss.branch[1] != 0) {
                        f = sel.node;
                        sel.val = true;
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

            for (size_t i = 0; i < path.size(); ++i) {
                if (path[i].val) {
                    itemset.insert(path[i].node.row());
                }
            }
        }
    };

    /**
     * Returns an iterator to the first instance,
     * which is viewed as a collection of item numbers.
     * Supports binary ZDDs only.
     * @return iterator to the first instance.
     */
    const_iterator begin() const {
        return const_iterator(*this, true);
    }

    /**
     * Returns an iterator to the element following the last instance.
     * Supports binary ZDDs only.
     * @return iterator to the instance following the last instance.
     */
    const_iterator end() const {
        return const_iterator(*this, false);
    }

    /**
     * Implements DdSpec.
     */
    int getRoot(nodeid& f) const {
        f = root_;
        return (f == 1) ? -1 : f.row();
    }

    /**
     * Implements DdSpec.
     */
    int getChild(nodeid& f, int level, int value) const {
        assert(level > 0 && level == f.row());
        assert(0 <= value && value < 2);
        f = child(f, value);
        return (f.row() > 0) ? f.row() : -f.col();
    }

    /**
     * Implements DdSpec.
     */
    size_t hashCode(nodeid const& f) const {
        return f.hash();
    }

    /**
     * Dumps the node table in Sapporo ZDD format.
     * Works only for binary DDs.
     * @param os the output stream.
     */
    void dumpSapporo(std::ostream& os) const {
        int const n = diagram->numRows() - 1;
        size_t const l = size();

        os << "_i " << n << "\n";
        os << "_o 1\n";
        os << "_n " << l << "\n";

        tdzdd::DataTable<size_t> nodeId(diagram->numRows());
        size_t k = 0;

        for (int i = 1; i <= n; ++i) {
            size_t const m = (*diagram)[i].size();
            Node<double> const* p = (*diagram)[i].data();
            nodeId[i].resize(m);

            for (size_t j = 0; j < m; ++j, ++p) {
                k += 2;
                nodeId[i][j] = k;
                os << k << " " << i;

                for (int c = 0; c <= 1; ++c) {
                    nodeid fc = p->branch[c];
                    if (fc == 0) {
                        os << " F";
                    }
                    else if (fc == 1) {
                        os << " T";
                    }
                    else {
                        os << " " << nodeId[fc.row()][fc.col()];
                    }
                }

                os << "\n";
            }

            tdzdd::MyVector<int> const& levels = diagram->lowerLevels(i);
            for (int const* t = levels.begin(); t != levels.end(); ++t) {
                nodeId[*t].clear();
            }
        }

        os << nodeId[root_.row()][root_.col()] << "\n";
        assert(k == l * 2);
    }
};
