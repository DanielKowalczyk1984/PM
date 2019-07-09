#ifndef NODE_BDD_REDUCER_HPP
#define NODE_BDD_REDUCER_HPP

#include <cassert>
#include <cmath>
#include <ostream>
#include <stdexcept>

#include <node_duration.hpp>
#include <NodeBddTable.hpp>
#include "util/MyHashTable.hpp"
// #include "tdzdd/util/MyList.hpp"
#include "util/MyVector.hpp"

template<typename T, bool BDD, bool ZDD>
class DdReducer {
    NodeTableEntity<T>& input;
    TableHandler<T> oldDiagram;
    TableHandler<T> newDiagram;
    NodeTableEntity<T>& output;
    MyVector<MyVector<NodeId> > newIdTable;
    MyVector<MyVector<NodeId*> > rootPtr;
    int counter = 1;

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
    explicit DdReducer(TableHandler<T>& diagram, bool useMP = false) :
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
            T* const tt = input[i].data();

            for (size_t j = 0; j < m; ++j) {
                for (int b = 0; b < 2; ++b) {
                    NodeId& f = tt[j].branch[b];
                    if (f.row() == 0) continue;

                    NodeId f0 = input.child(f, 0);
                    NodeId deletable = BDD ? f0 : 0;
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
    void setRoot(NodeId& root) {
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
        T* const tt = input[i].data();

        MyVector<NodeId>& newId = newIdTable[i];
        newId.resize(m);

        for (size_t j = m - 1; j + 1 > 0; --j) {
            NodeId& f0 = tt[j].branch[0];
            NodeId& f1 = tt[j].branch[1];

            if (f0.row() != 0) f0 = newIdTable[f0.row()][f0.col()];
            if (f1.row() != 0) f1 = newIdTable[f1.row()][f1.col()];

            if ((BDD && f1 == f0) || (ZDD && f1 == 0)) {
                newId[j] = f0;
            } else {
                newId[j] = NodeId(counter + 1, m); // tail of f0-equivalent list
            }
        }

        {
            MyVector<int> const& levels = input.lowerLevels(counter);
            for (int const* t = levels.begin(); t != levels.end(); ++t) {
                newIdTable[*t].clear();
            }
        }
        size_t mm = 0;

        for (size_t j = 0; j < m; ++j) {
            assert(newId[j].row() <= counter + 1);
            if (newId[j].row() <= counter){ continue;}

            NodeId& g0 = tt[j].branch[0];
            assert(newId[j].row() == counter + 1);
            newId[j] = NodeId(counter, mm++, g0.hasEmpty());
        }

        if (!BDD) {
            MyVector<int> const& levels = input.lowerLevels(counter);
            for (int const* t = levels.begin(); t != levels.end(); ++t) {
                input[*t].clear();
            }
        }

        if(mm > 0u) {
            output.initRow(counter, mm);
            T* nt = output[counter].data();

            for (size_t j = 0; j < m; ++j) {
                NodeId const& f0 = tt[j].branch[0];
                NodeId const& f1 = tt[j].branch[1];

                if ((BDD && f1 == f0) || (ZDD && f1 == 0)) { // forwarded
                    assert(newId[j].row() < counter);
                } else {
                    assert(newId[j].row() == counter);
                    size_t k = newId[j].col();
                    nt[k] = tt[j];
                    nt[k].set_head_node();
                    nt[k].child[0] = &(output.node(nt[k].branch[0]));
                    nt[k].child[1] = &(output.node(nt[k].branch[1]));
                }
            }

            counter++;
        }

        for (size_t k = 0; k < rootPtr[i].size(); ++k) {
            NodeId& root = *rootPtr[i][k];
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
            MyHashTable<Node<T> const*> uniq(m * 2);

            for (size_t j = 0; j < m; ++j) {
                Node<T>* const p0 = input[i].data();
                Node<T>& f = input[i][j];

                // make f canonical
                NodeId& f0 = f.branch[0];
                f0 = newIdTable[f0.row()][f0.col()];
                NodeId deletable = BDD ? f0 : 0;
                bool del = BDD || ZDD || (f0 == 0);
                for (int b = 1; b < 2; ++b) {
                    NodeId& ff = f.branch[b];
                    ff = newIdTable[ff.row()][ff.col()];
                    if (ff != deletable) del = false;
                }

                if (del) { // f is redundant
                    newIdTable[i][j] = f0;
                }
                else {
                    Node<T> const* pp = uniq.add(&f);

                    if (pp == &f) {
                        newIdTable[i][j] = NodeId(i, jj++, f0.hasEmpty());
                    }
                    else {
                        newIdTable[i][j] = newIdTable[i][pp - p0];
                    }
                }
            }
        }

        MyVector<int> const& levels = input.lowerLevels(i);
        for (int const* t = levels.begin(); t != levels.end(); ++t) {
            newIdTable[*t].clear();
        }

        output.initRow(i, jj);

        for (size_t j = 0; j < m; ++j) {
            NodeId const& ff = newIdTable[i][j];
            if (ff.row() == i) {
                output[i][ff.col()] = input[i][j];
                output[i][ff.col()].set_head_node();
                output[i][ff.col()].child[0] = &(output.node(output[i][ff.col()].branch[0]));
                output[i][ff.col()].child[1] = &(output.node(output[i][ff.col()].branch[1]));
            }
        }

        input[i].clear();

        for (size_t k = 0; k < rootPtr[i].size(); ++k) {
            NodeId& root = *rootPtr[i][k];
            root = newIdTable[i][root.col()];
        }
    }
};




#endif // NODE_BDD_REDUCER_HPP
