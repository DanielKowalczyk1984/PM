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

#ifndef NODE_BDD_BUILDER_HPP
#define NODE_BDD_BUILDER_HPP

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <unordered_set>

#include "NodeBdd.hpp"
#include "NodeBddSpec.hpp"
#include "NodeBddSweeper.hpp"
#include "NodeBddTable.hpp"
#include "NodeBranchId.hpp"
#include "util/MemoryPool.hpp"
#include "util/MyHashTable.hpp"
#include "util/MyList.hpp"
#include "util/MyVector.hpp"

class BuilderBase {
   protected:
    static int const headerSize = 1;

    /* SpecNode
     * ┌────────┬────────┬────────┬─────
     * │ srcPtr │state[0]│state[1]│ ...
     * │ nodeId │        │        │
     * └────────┴────────┴────────┴─────
     */
    union SpecNode {
        NodeId* srcPtr;
        int64_t code;
    };

    static NodeId*& srcPtr(SpecNode* p) { return p->srcPtr; }

    static int64_t& code(SpecNode* p) { return p->code; }

    static NodeId& nodeId(SpecNode* p) {
        return *reinterpret_cast<NodeId*>(&p->code);
    }

    static void* state(SpecNode* p) { return p + headerSize; }

    static void const* state(SpecNode const* p) { return p + headerSize; }

    static size_t getSpecNodeSize(int n) {
        if (n < 0) {
            throw std::runtime_error("storage size is not initialized!!!");
        }
        return headerSize + (n + sizeof(SpecNode) - 1) / sizeof(SpecNode);
    }

    template <typename SPEC>
    struct Hasher {
        SPEC const& spec;
        int const   level;

        Hasher(SPEC const& _spec, int _level) : spec(_spec), level(_level) {}

        size_t operator()(SpecNode const* p) const {
            return spec.hash_code(state(p), level);
        }

        size_t operator()(SpecNode const* p, SpecNode const* q) const {
            return spec.equal_to(state(p), state(q), level);
        }
    };
};

/**
 * Basic breadth-first DD builder.
 */
template <typename S, typename T = NodeBdd<double>>
class DdBuilder : BuilderBase {
    using Spec = S;
    using UniqTable = std::unordered_set<SpecNode*, Hasher<Spec>, Hasher<Spec>>;
    static int const AR = Spec::ARITY;

    Spec                spec;
    int const           specNodeSize;
    NodeTableEntity<T>& output;
    DdSweeper<T>        sweeper;

    std::vector<MyList<SpecNode>> spec_node_table;

    std::vector<char>         oneStorage;
    void* const               one;
    std::vector<NodeBranchId> oneSrcPtr;

    void init(int n) {
        spec_node_table.resize(n + 1);
        if (n >= output.numRows()) {
            output.setNumRows(n + 1);
        }
        oneSrcPtr.clear();
    }

   public:
    DdBuilder(Spec const& _spec, TableHandler<T>& _output, int n = 0)
        : spec(_spec),
          specNodeSize(getSpecNodeSize(_spec.datasize())),
          output(*_output),
          sweeper(this->output, oneSrcPtr),
          oneStorage(_spec.datasize()),
          one(oneStorage.data()) {
        if (n >= 1) {
            init(n);
        }
    }

    ~DdBuilder() {
        if (!oneSrcPtr.empty()) {
            spec.destruct(one);
            oneSrcPtr.clear();
        }
    }

    DdBuilder<S, T>(const DdBuilder<S, T>&) = default;
    DdBuilder<S, T>& operator=(const DdBuilder<S, T>&) = default;
    DdBuilder<S, T>(DdBuilder<S, T>&&) noexcept = default;
    DdBuilder<S, T>& operator=(DdBuilder<S, T>&&) noexcept = default;

    /**
     * Schedules a top-down event.
     * @param fp result storage.
     * @param level node level of the event.
     * @param s node state of the event.
     */
    void schedule(NodeId* fp, int level, void* s) {
        SpecNode* p0 = spec_node_table[level].alloc_front(specNodeSize);
        spec.get_copy(state(p0), s);
        srcPtr(p0) = fp;
    }

    /**
     * Initializes the builder.
     * @param root result storage.
     */
    int initialize(NodeId& root) {
        sweeper.setRoot(root);
        std::vector<char> tmp(spec.datasize());
        void* const       tmpState = tmp.data();
        int               n = spec.get_root(tmpState);

        if (n <= 0) {
            root = n ? 1 : 0;
            n = 0;
        } else {
            init(n);
            schedule(&root, n, tmpState);
        }

        spec.destruct(tmpState);
        if (!oneSrcPtr.empty()) {
            spec.destruct(one);
            oneSrcPtr.clear();
        }
        return n;
    }

    /**
     * Builds one level.
     * @param i level.
     */
    void construct(int i) {
        assert(0 < i && size_t(i) < spec_node_table.size());

        MyList<SpecNode>& spec_nodes = spec_node_table[i];
        size_t            j0 = output[i].size();
        size_t            m = j0;
        int               lowestChild = i - 1;
        size_t            deadCount = 0;

        {
            Hasher<Spec> hasher(spec, i);
            UniqTable    uniq(spec_nodes.size() * 2, hasher, hasher);

            for (auto p : spec_nodes) {
                // SpecNode*& p0 = uniq.add(p);
                auto aux = uniq.insert(p);

                if (aux.second) {
                    nodeId(p) = *srcPtr(p) = NodeId(i, m++);
                } else {
                    auto p0 = *(aux.first);
                    switch (spec.merge_states(state(p0), state(p))) {
                        case 1:
                            nodeId(p0) = 0;  // forward to 0-terminal
                            nodeId(p) = *srcPtr(p) = NodeId(i, m++);
                            p0 = p;
                            break;
                        case 2:
                            *srcPtr(p) = 0;
                            nodeId(p) = 1;  // unused
                            break;
                        default:
                            *srcPtr(p) = nodeId(p0);
                            nodeId(p) = 1;  // unused
                            break;
                    }
                }
            }
            //#ifdef DEBUG
            //            MessageHandler mh;
            //            mh << "table_size[" << i << "] = " << uniq.tableSize()
            //            << "\n";
            //#endif
        }

        output[i].resize(m);
        T* const  output_data = output[i].data();
        size_t    jj = j0;
        SpecNode* pp = spec_node_table[i - 1].alloc_front(specNodeSize);

        for (; !spec_nodes.empty(); spec_nodes.pop_front()) {
            SpecNode* p = spec_nodes.front();
            T&        q = output_data[jj];

            if (nodeId(p) == 1) {
                spec.destruct(state(p));
                continue;
            }

            bool allZero = true;

            for (int b = 0; b < AR; ++b) {
                if (nodeId(p) == 0) {
                    q[b] = 0;
                    continue;
                }

                spec.get_copy(state(pp), state(p));
                int ii = spec.get_child(state(pp), i, b);

                if (ii == 0) {
                    q[b] = 0;
                    spec.destruct(state(pp));
                } else if (ii < 0) {
                    if (oneSrcPtr.empty()) {  // the first 1-terminal candidate
                        spec.get_copy(one, state(pp));
                        q[b] = 1;
                        oneSrcPtr.emplace_back(i, jj, b);
                    } else {
                        switch (spec.merge_states(one, state(pp))) {
                            case 1:
                                while (!oneSrcPtr.empty()) {
                                    NodeBranchId const& nbi = oneSrcPtr.back();
                                    assert(nbi.row >= i);
                                    output[nbi.row][nbi.col][nbi.val] = 0;
                                    oneSrcPtr.pop_back();
                                }
                                spec.destruct(one);
                                spec.get_copy(one, state(pp));
                                q[b] = 1;
                                oneSrcPtr.emplace_back(i, jj, b);
                                break;
                            case 2:
                                q[b] = 0;
                                break;
                            default:
                                q[b] = 1;
                                oneSrcPtr.emplace_back(i, jj, b);
                                break;
                        }
                    }
                    spec.destruct(state(pp));
                    allZero = false;
                } else if (ii == i - 1) {
                    srcPtr(pp) = &q[b];
                    pp = spec_node_table[ii].alloc_front(specNodeSize);
                    allZero = false;
                } else {
                    assert(ii < i - 1);
                    SpecNode* ppp =
                        spec_node_table[ii].alloc_front(specNodeSize);
                    spec.get_copy(state(ppp), state(pp));
                    spec.destruct(state(pp));
                    srcPtr(ppp) = &q[b];
                    if (ii < lowestChild) {
                        lowestChild = ii;
                    }
                    allZero = false;
                }
            }

            spec.destruct(state(p));
            ++jj;
            if (allZero) {
                ++deadCount;
            }
        }

        spec_node_table[i - 1].pop_front();
        // spec.destructLevel(i);
        sweeper.update(i, lowestChild, deadCount);
    }
};

/**
 * Multi-threaded breadth-first DD builder.
 */

/**
 * Breadth-first ZDD subset builder.
 */
template <typename T, typename S>
class ZddSubsetter : BuilderBase {
    // typedef typename std::remove_const<typename
    // std::remove_reference<S>::type>::type Spec;
    using Spec = S;
    using UniqTable = std::unordered_set<SpecNode*, Hasher<Spec>, Hasher<Spec>>;
    static int const AR = Spec::ARITY;

    Spec                              spec;
    int const                         specNodeSize;
    NodeTableEntity<T> const&         input;
    NodeTableEntity<T>&               output;
    DataTable<MyListOnPool<SpecNode>> work;
    DdSweeper<T>                      sweeper;

    std::vector<char>         oneStorage;
    void* const               one;
    std::vector<NodeBranchId> oneSrcPtr;

    MemoryPools pools;

   public:
    ZddSubsetter(TableHandler<T> const& _input,
                 Spec const&            s,
                 TableHandler<T>&       _output)
        : spec(s),
          specNodeSize(getSpecNodeSize(spec.datasize())),
          input(*_input),
          output(*_output),
          work(_input->numRows()),
          sweeper(this->output, oneSrcPtr),
          oneStorage(spec.datasize()),
          one(oneStorage.data()) {}

    ~ZddSubsetter() {
        if (!oneSrcPtr.empty()) {
            spec.destruct(one);
            oneSrcPtr.clear();
        }
    }

    ZddSubsetter<T, S>(const ZddSubsetter<T, S>&) = default;
    ZddSubsetter<T, S>& operator=(const ZddSubsetter<T, S>&) = default;
    ZddSubsetter<T, S>(ZddSubsetter<T, S>&&) noexcept = default;
    ZddSubsetter<T, S>& operator=(ZddSubsetter<T, S>&&) noexcept = default;

    /**
     * Initializes the builder.
     * @param root the root node.
     */
    int initialize(NodeId& root) {
        sweeper.setRoot(root);
        std::vector<char> tmp(spec.datasize());
        void* const       tmpState = tmp.data();
        int               n = spec.get_root(tmpState);

        int k = (root == 1) ? -1 : root.row();

        while (n != 0 && k != 0 && n != k) {
            if (n < k) {
                assert(k >= 1);
                k = downTable(root, 0, n);
            } else {
                assert(n >= 1);
                n = downSpec(tmpState, n, 0, k);
            }
        }

        if (n <= 0 || k <= 0) {
            assert(n == 0 || k == 0 || (n == -1 && k == -1));
            root = NodeId(0, n != 0 && k != 0);
            n = 0;
        } else {
            assert(n == k);
            assert(n == root.row());

            pools.resize(n + 1);
            work[n].resize(input[n].size());

            SpecNode* p0 =
                work[n][root.col()].alloc_front(pools[n], specNodeSize);
            spec.get_copy(state(p0), tmpState);
            srcPtr(p0) = &root;
        }

        spec.destruct(tmpState);
        output.init(n + 1);
        if (!oneSrcPtr.empty()) {
            spec.destruct(one);
            oneSrcPtr.clear();
        }
        return n;
    }

    /**
     * Builds one level.
     * @param i level.
     */
    void subset(int i) {
        assert(0 < i && i < output.numRows());
        assert(output.numRows() - pools.size() == 0);

        Hasher<Spec> const hasher(spec, i);
        std::vector<char>  tmp(spec.datasize());
        void* const        tmpState = tmp.data();
        size_t const       m = input[i].size();
        size_t             mm = 0;
        int                lowestChild = i - 1;
        size_t             deadCount = 0;

        if (work[i].empty()) {
            work[i].resize(m);
        }
        assert(work[i].size() == m);

        for (size_t j = 0; j < m; ++j) {
            MyListOnPool<SpecNode>& list = work[i][j];
            size_t                  n = list.size();

            if (n >= 2) {
                UniqTable uniq(n * 2, hasher, hasher);

                for (auto p : list) {
                    auto aux = uniq.insert(p);
                    // SpecNode*& p0 = uniq.add(p);

                    if (aux.second) {
                        nodeId(p) = *srcPtr(p) = NodeId(i, mm++);
                    } else {
                        auto p0 = *(aux.first);
                        switch (spec.merge_states(state(p0), state(p))) {
                            case 1:
                                nodeId(p0) = 0;  // forward to 0-terminal
                                nodeId(p) = *srcPtr(p) = NodeId(i, mm++);
                                p0 = p;
                                break;
                            case 2:
                                *srcPtr(p) = 0;
                                nodeId(p) = 1;  // unused
                                break;
                            default:
                                *srcPtr(p) = nodeId(p0);
                                nodeId(p) = 1;  // unused
                                break;
                        }
                    }
                }
            } else if (n == 1) {
                SpecNode* p = list.front();
                nodeId(p) = *srcPtr(p) = NodeId(i, mm++);
            }
        }

        output.initRow(i, mm);
        T* const output_data = output[i].data();
        size_t   jj = 0;

        for (size_t j = 0; j < m; ++j) {
            MyListOnPool<SpecNode>& list = work[i][j];

            for (auto p : list) {
                T& q = output_data[jj];

                if (nodeId(p) == 1) {
                    spec.destruct(state(p));
                    continue;
                }

                bool allZero = true;

                for (int b = 0; b < AR; ++b) {
                    if (nodeId(p) == 0) {
                        q[b] = 0;
                        continue;
                    }

                    NodeId f(i, j);
                    spec.get_copy(tmpState, state(p));
                    int kk = downTable(f, b, i - 1);
                    int ii = downSpec(tmpState, i, b, kk);

                    while (ii != 0 && kk != 0 && ii != kk) {
                        if (ii < kk) {
                            assert(kk >= 1);
                            kk = downTable(f, 0, ii);
                        } else {
                            assert(ii >= 1);
                            ii = downSpec(tmpState, ii, 0, kk);
                        }
                    }

                    if (ii <= 0 || kk <= 0) {
                        if (ii == 0 || kk == 0) {
                            q[b] = 0;
                        } else {
                            if (oneSrcPtr.empty()) {  // the first 1-terminal
                                                      // candidate
                                spec.get_copy(one, tmpState);
                                q[b] = 1;
                                oneSrcPtr.emplace_back(i, jj, b);
                            } else {
                                switch (spec.merge_states(one, tmpState)) {
                                    case 1:
                                        while (!oneSrcPtr.empty()) {
                                            NodeBranchId const& nbi =
                                                oneSrcPtr.back();
                                            assert(nbi.row >= i);
                                            output[nbi.row][nbi.col][nbi.val] =
                                                0;
                                            oneSrcPtr.pop_back();
                                        }
                                        spec.destruct(one);
                                        spec.get_copy(one, tmpState);
                                        q[b] = 1;
                                        oneSrcPtr.emplace_back(i, jj, b);
                                        break;
                                    case 2:
                                        q[b] = 0;
                                        break;
                                    default:
                                        q[b] = 1;
                                        oneSrcPtr.emplace_back(i, jj, b);
                                        break;
                                }
                            }
                            allZero = false;
                        }
                    } else {
                        assert(ii == f.row() && ii == kk && ii < i);
                        if (work[ii].empty()) {
                            work[ii].resize(input[ii].size());
                        }
                        SpecNode* pp = work[ii][f.col()].alloc_front(
                            pools[ii], specNodeSize);
                        spec.get_copy(state(pp), tmpState);
                        srcPtr(pp) = &q[b];
                        if (ii < lowestChild) {
                            lowestChild = ii;
                        }
                        allZero = false;
                    }

                    spec.destruct(tmpState);
                }

                spec.destruct(state(p));
                ++jj;
                if (allZero) {
                    ++deadCount;
                }
            }
        }

        work[i].clear();
        pools[i].clear();
        // spec.destructLevel(i);
        sweeper.update(i, lowestChild, deadCount);
    }

   private:
    int downTable(NodeId& f, int b, int zerosupLevel) const {
        if (zerosupLevel < 0) {
            zerosupLevel = 0;
        }

        f = input.child(f, b);
        while (f.row() > zerosupLevel) {
            f = input.child(f, 0);
        }
        return (f == 1) ? -1 : f.row();
    }

    int downSpec(void* p, int level, int b, int zerosupLevel) {
        if (zerosupLevel < 0) {
            zerosupLevel = 0;
        }
        assert(level > zerosupLevel);

        int i = spec.get_child(p, level, b);
        while (i > zerosupLevel) {
            i = spec.get_child(p, i, 0);
        }
        return i;
    }
};

#endif  // NODE_BDD_BUILDER_HPP
