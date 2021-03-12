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

#ifndef NODE_BDD_SWEEPER_HPP
#define NODE_BDD_SWEEPER_HPP

#include <cassert>
#include <ostream>
#include <vector>

#include "NodeBdd.hpp"
#include "NodeBddTable.hpp"
#include "NodeBranchId.hpp"
#include "util/MyVector.hpp"

/**
 * On-the-fly DD cleaner.
 * Removes the nodes that are identified as equivalent to the 0-terminal
 * while top-down DD construction.
 */
template <typename T = NodeBdd<double>>
class DdSweeper {
    static size_t const SWEEP_RATIO = 20;

    NodeTableEntity<T>&        diagram;
    std::vector<NodeBranchId>* oneSrcPtr;

    std::vector<int>    sweepLevel;
    std::vector<size_t> deadCount;
    size_t              allCount;
    size_t              maxCount;
    NodeId*             rootPtr;

   public:
    /**
     * Constructor.
     * @param diagram the diagram to sweep.
     */
    explicit DdSweeper(NodeTableEntity<T>& _diagram)
        : diagram(_diagram),
          oneSrcPtr(nullptr),
          allCount(0),
          maxCount(0),
          rootPtr(nullptr) {}

    /**
     * Constructor.
     * @param diagram the diagram to sweep.
     * @param oneSrcPtr collection of node branch IDs.
     */
    DdSweeper(NodeTableEntity<T>&        _diagram,
              std::vector<NodeBranchId>& _oneSrcPtr)
        : diagram(_diagram),
          oneSrcPtr(&_oneSrcPtr),
          allCount(0),
          maxCount(0),
          rootPtr(nullptr) {}

    /**
     * Set the root pointer.
     * @param root reference to the root ID storage.
     */
    void setRoot(NodeId& root) { rootPtr = &root; }

    /**
     * Updates status and sweeps the DD if necessary.
     * @param current current level.
     * @param child the level at which edges from this level are completed.
     * @param count the number of dead nodes at this level.
     */
    void update(int current, int child, size_t count) {
        assert(1 <= current);
        assert(0 <= child);
        if (current <= 1) {
            return;
        }

        if (size_t(current) >= sweepLevel.size()) {
            sweepLevel.resize(current + 1);
            deadCount.resize(current + 2);
        }

        for (int i = child; i <= current; ++i) {
            if (sweepLevel[i] > 0) {
                break;
            }
            sweepLevel[i] = current + 1;
        }

        deadCount[current] = count;
        allCount += diagram[current].size();

        int k = sweepLevel[current - 1];
        for (int i = sweepLevel[current]; i > k; --i) {
            deadCount[k] += deadCount[i];
            deadCount[i] = 0;
        }
        if (maxCount < allCount) {
            maxCount = allCount;
        }
        if (deadCount[k] * SWEEP_RATIO < maxCount) {
            return;
        }

        std::vector<std::vector<NodeId>> newId(diagram.numRows());

        for (int i = k; i < diagram.numRows(); ++i) {
            size_t m = diagram[i].size();
            newId[i].resize(m);

            size_t jj = 0;

            for (size_t j = 0; j < m; ++j) {
                T&   p = diagram[i][j];
                bool dead = true;

                for (int b = 0; b < 2; ++b) {
                    NodeId& f = p[b];
                    if (f.row() >= k) {
                        f = newId[f.row()][f.col()];
                    }
                    if (f != 0) {
                        dead = false;
                    }
                }

                if (dead) {
                    newId[i][j] = 0;
                } else {
                    newId[i][j] = NodeId(i, jj);
                    diagram[i][jj] = p;
                    ++jj;
                }
            }

            diagram[i].resize(jj);
        }

        if (oneSrcPtr) {
            for (auto& nbi : *oneSrcPtr) {
                if (nbi.row >= k) {
                    NodeId f = newId[nbi.row][nbi.col];
                    nbi.row = f.row();
                    nbi.col = f.col();
                }
            }
        }

        *rootPtr = newId[rootPtr->row()][rootPtr->col()];
        deadCount[k] = 0;
        allCount = diagram.size();
    }
};

#endif  // NODE_BDD_SWEEPER_HPP
