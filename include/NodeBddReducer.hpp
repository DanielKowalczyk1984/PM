#ifndef NODE_BDD_REDUCER_HPP
#define NODE_BDD_REDUCER_HPP

#include <cassert>
#include <cmath>
#include <ostream>
#include <stdexcept>

#include <node_duration.hpp>
#include <NodeBddTable.hpp>
#include "tdzdd/util/MyHashTable.hpp"
#include "tdzdd/util/MyList.hpp"
#include "tdzdd/util/MyVector.hpp"

template<typename T, bool BDD, bool ZDD>
class DdReducer {
    NodeTableEntity<T>& input;
    NodeTableHandler<T> oldDiagram;
    NodeTableHandler<T> newDiagram;
    NodeTableEntity<T>& output;
    tdzdd::MyVector<tdzdd::MyVector<nodeid> > newIdTable;
    tdzdd::MyVector<tdzdd::MyVector<nodeid*> > rootPtr;

    struct ReducNodeInfo {
        Node<T> children;
        size_t column;

        size_t hash() const {
            return children.hash();
        }

        bool operator==(ReducNodeInfo const& o) const {
            return children == o.children;
        }

        friend std::ostream& operator<<(std::ostream& os,
                                        ReducNodeInfo const& o) {
            return os << "(" << o.children << " -> " << o.column << ")";
        }
    };

    bool readyForSequentialReduction;

public:
    explicit DdReducer(NodeTableHandler<T>& diagram, bool useMP = false) :
            input(diagram.privateEntity()),
            oldDiagram(diagram),
            newDiagram(input.numRows()),
            output(newDiagram.privateEntity()),
            newIdTable(input.numRows()),
            rootPtr(input.numRows()),
            readyForSequentialReduction(false) {
        diagram = newDiagram;

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
        if (readyForSequentialReduction) return;

        for (int i = 2; i < input.numRows(); ++i) {
            size_t const m = input[i].size();
            Node<T>* const tt = input[i].data();

            for (size_t j = 0; j < m; ++j) {
                for (int b = 0; b < 2; ++b) {
                    nodeid& f = tt[j].branch[b];
                    if (f.row() == 0) continue;

                    nodeid f0 = input.child(f, 0);
                    nodeid deletable = BDD ? f0 : 0;
                    bool del = true;

                    for (int bb = (BDD || ZDD) ? 1 : 0; bb < 2; ++bb) {
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
    void setRoot(nodeid& root) {
        rootPtr[root.row()].push_back(&root);
    }

    /**
     * Reduces one level.
     * @param i level.
     * @param useMP use an algorithm for multiple processors.
     */
    void reduce(int i) {
        // if (2 == 2) {
        algorithmR(i);
        // }
        // else {
            // reduce_(i);
        // }
    }

private:
    /**
     * Reduces one level using Algorithm-R.
     * @param i level.
     */
    void algorithmR(int i) {
        makeReadyForSequentialReduction();
        size_t const m = input[i].size();
        Node<T>* const tt = input[i].data();
        nodeid const mark(i, m);

        tdzdd::MyVector<nodeid>& newId = newIdTable[i];
        newId.resize(m);

        for (size_t j = m - 1; j + 1 > 0; --j) {
            nodeid& f0 = tt[j].branch[0];
            nodeid& f1 = tt[j].branch[1];

            if (f0.row() != 0) f0 = newIdTable[f0.row()][f0.col()];
            if (f1.row() != 0) f1 = newIdTable[f1.row()][f1.col()];

            if ((BDD && f1 == f0) || (ZDD && f1 == 0)) {
                newId[j] = f0;
            }
            else {
                nodeid& f00 = input.child(f0, 0);
                nodeid& f01 = input.child(f0, 1);

                // if (f01 != mark) {        // the first touch from this level
                    f01 = mark;        // mark f0 as touched
                    newId[j] = nodeid(i + 1, m); // tail of f0-equivalent list
                // }
                // else {
                //     newId[j] = f00;         // next of f0-equivalent list
                // }
                f00 = nodeid(i + 1, j);  // new head of f0-equivalent list
            }
        }

        {
            tdzdd::MyVector<int> const& levels = input.lowerLevels(i);
            for (int const* t = levels.begin(); t != levels.end(); ++t) {
                newIdTable[*t].clear();
            }
        }
        size_t mm = 0;

        for (size_t j = 0; j < m; ++j) {
            nodeid const f(i, j);
            assert(newId[j].row() <= i + 1);
            if (newId[j].row() <= i){ continue;}

            for (size_t k = j; k < m;) { // for each g in f0-equivalent list
                assert(j <= k);
                nodeid const g(i, k);
                nodeid& g0 = tt[k].branch[0];
                nodeid& g1 = tt[k].branch[1];
                nodeid& g10 = input.child(g1, 0);
                nodeid& g11 = input.child(g1, 1);
                assert(g1 != mark);
                assert(newId[k].row() == i + 1);
                size_t next = newId[k].col();

                if (g11 != f) { // the first touch to g1 in f0-equivalent list
                    g11 = f; // mark g1 as touched
                    g10 = g; // record g as a canonical node for <f0,g1>
                    newId[k] = nodeid(i, mm++, g0.hasEmpty());
                }
                else {
                    g0 = g10;       // make a forward link
                    g1 = mark;      // mark g as forwarded
                    newId[k] = 0;
                }

                k = next;
            }
        }

        if (!BDD) {
            tdzdd::MyVector<int> const& levels = input.lowerLevels(i);
            for (int const* t = levels.begin(); t != levels.end(); ++t) {
                input[*t].clear();
            }
        }

        output.initRow(i, mm);
        Node<T>* nt = output[i].data();

        for (size_t j = 0; j < m; ++j) {
            nodeid const& f0 = tt[j].branch[0];
            nodeid const& f1 = tt[j].branch[1];

            if (f1 == mark) { // forwarded
                assert(f0.row() == i);
                assert(newId[j] == 0);
                newId[j] = newId[f0.col()];
            }
            else if ((BDD && f1 == f0) || (ZDD && f1 == 0)) { // forwarded
                assert(newId[j].row() < i);
            }
            else {
                assert(newId[j].row() == i);
                size_t k = newId[j].col();
                nt[k] = tt[j];
                nt[k].set_head_node();
                nt[k].child[0] = &(output.node(nt[k].branch[0]));
                nt[k].child[1] = &(output.node(nt[k].branch[1]));
            }
        }

        for (size_t k = 0; k < rootPtr[i].size(); ++k) {
            nodeid& root = *rootPtr[i][k];
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
            //MyList<ReducNodeInfo> rni;
            //MyHashTable<ReducNodeInfo const*> uniq(m * 2);
            tdzdd::MyHashTable<Node<T> const*> uniq(m * 2);

            for (size_t j = 0; j < m; ++j) {
                Node<T>* const p0 = input[i].data();
                Node<T>& f = input[i][j];

                // make f canonical
                nodeid& f0 = f.branch[0];
                f0 = newIdTable[f0.row()][f0.col()];
                nodeid deletable = BDD ? f0 : 0;
                bool del = BDD || ZDD || (f0 == 0);
                for (int b = 1; b < 2; ++b) {
                    nodeid& ff = f.branch[b];
                    ff = newIdTable[ff.row()][ff.col()];
                    if (ff != deletable) del = false;
                }

                if (del) { // f is redundant
                    newIdTable[i][j] = f0;
                }
                else {
                    Node<T> const* pp = uniq.add(&f);

                    if (pp == &f) {
                        newIdTable[i][j] = nodeid(i, jj++, f0.hasEmpty());
                    }
                    else {
                        newIdTable[i][j] = newIdTable[i][pp - p0];
                    }
                }
            }
        }

        tdzdd::MyVector<int> const& levels = input.lowerLevels(i);
        for (int const* t = levels.begin(); t != levels.end(); ++t) {
            newIdTable[*t].clear();
        }

        output.initRow(i, jj);

        for (size_t j = 0; j < m; ++j) {
            nodeid const& ff = newIdTable[i][j];
            if (ff.row() == i) {
                output[i][ff.col()] = input[i][j];
                output[i][ff.col()].set_head_node();
                output[i][ff.col()].child[0] = &(output.node(output[i][ff.col()].branch[0]));
                output[i][ff.col()].child[1] = &(output.node(output[i][ff.col()].branch[1]));
            }
        }

        input[i].clear();

        for (size_t k = 0; k < rootPtr[i].size(); ++k) {
            nodeid& root = *rootPtr[i][k];
            root = newIdTable[i][root.col()];
        }
    }
};




#endif // NODE_BDD_REDUCER_HPP
