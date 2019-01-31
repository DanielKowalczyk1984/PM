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

#include <cassert>
#include <cmath>
#include <ostream>
#include <stdexcept>

#include "NodeBddSweeper.hpp"
#include "node_duration.hpp"
#include "NodeBddTable.hpp"
#include "tdzdd/DdSpec.hpp"
#include "tdzdd/util/MemoryPool.hpp"
#include "tdzdd/util/MyHashTable.hpp"
#include "tdzdd/util/MyList.hpp"
#include "tdzdd/util/MyVector.hpp"


class DdBuilderBase {
protected:
    static int const headerSize = 1;

    /* SpecNode
     * ┌────────┬────────┬────────┬─────
     * │ srcPtr │state[0]│state[1]│ ...
     * │ nodeId │        │        │
     * └────────┴────────┴────────┴─────
     */
    union SpecNode {
        nodeid* srcPtr;
        int64_t code;
    };

    static nodeid*& srcPtr(SpecNode* p) {
        return p[0].srcPtr;
    }

    static int64_t& code(SpecNode* p) {
        return p[0].code;
    }

    static nodeid& nodeId(SpecNode* p) {
        return *reinterpret_cast<nodeid*>(&p[0].code);
    }

    static void* state(SpecNode* p) {
        return p + headerSize;
    }

    static void const* state(SpecNode const* p) {
        return p + headerSize;
    }

    static int getSpecNodeSize(int n) {
        if (n < 0)
            throw std::runtime_error("storage size is not initialized!!!");
        return headerSize + (n + sizeof(SpecNode) - 1) / sizeof(SpecNode);
    }

    template<typename SPEC>
    struct Hasher {
        SPEC const& spec;
        int const level;

        Hasher(SPEC const& spec, int level) :
                spec(spec), level(level) {
        }

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
template< typename S, typename T = double>
class DdBuilder: DdBuilderBase {
    typedef S Spec;
    typedef tdzdd::MyHashTable<SpecNode*,Hasher<Spec>,Hasher<Spec> > UniqTable;
    static int const AR = Spec::ARITY;

    Spec spec;
    int const specNodeSize;
    NodeTableEntity<>& output;
    DdSweeper<T> sweeper;

    tdzdd::MyVector<tdzdd::MyList<SpecNode> > snodeTable;

    tdzdd::MyVector<char> oneStorage;
    void* const one;
    tdzdd::MyVector<tdzdd::NodeBranchId> oneSrcPtr;

    void init(int n) {
        snodeTable.resize(n + 1);
        if (n >= output.numRows()) output.setNumRows(n + 1);
        oneSrcPtr.clear();
    }

public:
    DdBuilder(Spec const& spec, NodeTableHandler<T>& output, int n = 0) :
            spec(spec),
            specNodeSize(getSpecNodeSize(spec.datasize())),
            output(output.privateEntity()),
            sweeper(this->output, oneSrcPtr),
            oneStorage(spec.datasize()),
            one(oneStorage.data()) {
        if (n >= 1) init(n);
    }

    ~DdBuilder() {
        if (!oneSrcPtr.empty()) {
            spec.destruct(one);
            oneSrcPtr.clear();
        }
    }

    /**
     * Schedules a top-down event.
     * @param fp result storage.
     * @param level node level of the event.
     * @param s node state of the event.
     */
    void schedule(nodeid* fp, int level, void* s) {
        SpecNode* p0 = snodeTable[level].alloc_front(specNodeSize);
        spec.get_copy(state(p0), s);
        srcPtr(p0) = fp;
    }

    /**
     * Initializes the builder.
     * @param root result storage.
     */
    int initialize(nodeid& root) {
        sweeper.setRoot(root);
        tdzdd::MyVector<char> tmp(spec.datasize());
        void* const tmpState = tmp.data();
        int n = spec.get_root(tmpState);

        if (n <= 0) {
            root = n ? 1 : 0;
            n = 0;
        }
        else {
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
        assert(0 < i && size_t(i) < snodeTable.size());

        tdzdd::MyList<SpecNode> &snodes = snodeTable[i];
        size_t j0 = output[i].size();
        size_t m = j0;
        int lowestChild = i - 1;
        size_t deadCount = 0;

        {
            Hasher<Spec> hasher(spec, i);
            UniqTable uniq(snodes.size() * 2, hasher, hasher);

            for (tdzdd::MyList<SpecNode>::iterator t = snodes.begin();
                    t != snodes.end(); ++t) {
                SpecNode* p = *t;
                SpecNode*& p0 = uniq.add(p);

                if (p0 == p) {
                    nodeId(p) = *srcPtr(p) = nodeid(i, m++);
                }
                else {
                    switch (spec.merge_states(state(p0), state(p))) {
                    case 1:
                        nodeId(p0) = 0; // forward to 0-terminal
                        nodeId(p) = *srcPtr(p) = nodeid(i, m++);
                        p0 = p;
                        break;
                    case 2:
                        *srcPtr(p) = 0;
                        nodeId(p) = 1; // unused
                        break;
                    default:
                        *srcPtr(p) = nodeId(p0);
                        nodeId(p) = 1; // unused
                        break;
                    }
                }
            }
//#ifdef DEBUG
//            MessageHandler mh;
//            mh << "table_size[" << i << "] = " << uniq.tableSize() << "\n";
//#endif
        }

        output[i].resize(m);
        Node<T>* const outi = output[i].data();
        size_t jj = j0;
        SpecNode* pp = snodeTable[i - 1].alloc_front(specNodeSize);

        for (; !snodes.empty(); snodes.pop_front()) {
            SpecNode* p = snodes.front();
            Node<T>& q = outi[jj];

            if (nodeId(p) == 1) {
                spec.destruct(state(p));
                continue;
            }

            bool allZero = true;

            for (int b = 0; b < AR; ++b) {
                if (nodeId(p) == 0) {
                    q.branch[b] = 0;
                    continue;
                }

                spec.get_copy(state(pp), state(p));
                int ii = spec.get_child(state(pp), i, b);

                if (ii == 0) {
                    q.branch[b] = 0;
                    spec.destruct(state(pp));
                }
                else if (ii < 0) {
                    if (oneSrcPtr.empty()) { // the first 1-terminal candidate
                        spec.get_copy(one, state(pp));
                        q.branch[b] = 1;
                        oneSrcPtr.push_back(tdzdd::NodeBranchId(i, jj, b));
                    }
                    else {
                        switch (spec.merge_states(one, state(pp))) {
                        case 1:
                            while (!oneSrcPtr.empty()) {
                                tdzdd::NodeBranchId const& nbi = oneSrcPtr.back();
                                assert(nbi.row >= i);
                                output[nbi.row][nbi.col].branch[nbi.val] = 0;
                                oneSrcPtr.pop_back();
                            }
                            spec.destruct(one);
                            spec.get_copy(one, state(pp));
                            q.branch[b] = 1;
                            oneSrcPtr.push_back(tdzdd::NodeBranchId(i, jj, b));
                            break;
                        case 2:
                            q.branch[b] = 0;
                            break;
                        default:
                            q.branch[b] = 1;
                            oneSrcPtr.push_back(tdzdd::NodeBranchId(i, jj, b));
                            break;
                        }
                    }
                    spec.destruct(state(pp));
                    allZero = false;
                }
                else if (ii == i - 1) {
                    srcPtr(pp) = &q.branch[b];
                    pp = snodeTable[ii].alloc_front(specNodeSize);
                    allZero = false;
                }
                else {
                    assert(ii < i - 1);
                    SpecNode* ppp = snodeTable[ii].alloc_front(specNodeSize);
                    spec.get_copy(state(ppp), state(pp));
                    spec.destruct(state(pp));
                    srcPtr(ppp) = &q.branch[b];
                    if (ii < lowestChild) lowestChild = ii;
                    allZero = false;
                }
            }

            spec.destruct(state(p));
            ++jj;
            if (allZero) ++deadCount;
        }

        snodeTable[i - 1].pop_front();
        spec.destructLevel(i);
        sweeper.update(i, lowestChild, deadCount);
    }
};

/**
 * Multi-threaded breadth-first DD builder.
 */

/**
 * Breadth-first ZDD subset builder.
 */
template<typename T, typename S>
class ZddSubsetter: DdBuilderBase {
//typedef typename std::remove_const<typename std::remove_reference<S>::type>::type Spec;
    typedef S Spec;
    typedef tdzdd::MyHashTable<SpecNode*,Hasher<Spec>,Hasher<Spec> > UniqTable;
    static int const AR = Spec::ARITY;

    Spec spec;
    int const specNodeSize;
    NodeTableEntity<T> const& input;
    NodeTableEntity<T>& output;
    tdzdd::DataTable<tdzdd::MyListOnPool<SpecNode> > work;
    tdzdd::DdSweeper<AR> sweeper;

    tdzdd::MyVector<char> oneStorage;
    void* const one;
    tdzdd::MyVector<tdzdd::NodeBranchId> oneSrcPtr;

    tdzdd::MemoryPools pools;

public:
    ZddSubsetter(NodeTableHandler<T> const& input, Spec const& s,
                 NodeTableHandler<T>& output) :
            spec(s),
            specNodeSize(getSpecNodeSize(spec.datasize())),
            input(*input),
            output(output.privateEntity()),
            work(input->numRows()),
            sweeper(this->output, oneSrcPtr),
            oneStorage(spec.datasize()),
            one(oneStorage.data()) {
    }

    ~ZddSubsetter() {
        if (!oneSrcPtr.empty()) {
            spec.destruct(one);
            oneSrcPtr.clear();
        }
    }

    /**
     * Initializes the builder.
     * @param root the root node.
     */
    int initialize(nodeid& root) {
        sweeper.setRoot(root);
        tdzdd::MyVector<char> tmp(spec.datasize());
        void* const tmpState = tmp.data();
        int n = spec.get_root(tmpState);

        int k = (root == 1) ? -1 : root.row();

        while (n != 0 && k != 0 && n != k) {
            if (n < k) {
                assert(k >= 1);
                k = downTable(root, 0, n);
            }
            else {
                assert(n >= 1);
                n = downSpec(tmpState, n, 0, k);
            }
        }

        if (n <= 0 || k <= 0) {
            assert(n == 0 || k == 0 || (n == -1 && k == -1));
            root = nodeid(0, n != 0 && k != 0);
            n = 0;
        }
        else {
            assert(n == k);
            assert(n == root.row());

            pools.resize(n + 1);
            work[n].resize(input[n].size());

            SpecNode* p0 = work[n][root.col()].alloc_front(pools[n],
                    specNodeSize);
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
        tdzdd::MyVector<char> tmp(spec.datasize());
        void* const tmpState = tmp.data();
        size_t const m = input[i].size();
        size_t mm = 0;
        int lowestChild = i - 1;
        size_t deadCount = 0;

        if (work[i].empty()) work[i].resize(m);
        assert(work[i].size() == m);

        for (size_t j = 0; j < m; ++j) {
            tdzdd::MyListOnPool<SpecNode> &list = work[i][j];
            size_t n = list.size();

            if (n >= 2) {
                UniqTable uniq(n * 2, hasher, hasher);

                for (tdzdd::MyListOnPool<SpecNode>::iterator t = list.begin();
                        t != list.end(); ++t) {
                    SpecNode* p = *t;
                    SpecNode*& p0 = uniq.add(p);

                    if (p0 == p) {
                        nodeId(p) = *srcPtr(p) = nodeid(i, mm++);
                    }
                    else {
                        switch (spec.merge_states(state(p0), state(p))) {
                        case 1:
                            nodeId(p0) = 0; // forward to 0-terminal
                            nodeId(p) = *srcPtr(p) = nodeid(i, mm++);
                            p0 = p;
                            break;
                        case 2:
                            *srcPtr(p) = 0;
                            nodeId(p) = 1; // unused
                            break;
                        default:
                            *srcPtr(p) = nodeId(p0);
                            nodeId(p) = 1; // unused
                            break;
                        }
                    }
                }
            }
            else if (n == 1) {
                SpecNode* p = list.front();
                nodeId(p) = *srcPtr(p) = nodeid(i, mm++);
            }
        }

        output.initRow(i, mm);
        Node<T>* const outi = output[i].data();
        size_t jj = 0;

        for (size_t j = 0; j < m; ++j) {
            tdzdd::MyListOnPool<SpecNode> &list = work[i][j];

            for (tdzdd::MyListOnPool<SpecNode>::iterator t = list.begin();
                    t != list.end(); ++t) {
                SpecNode* p = *t;
                Node<T>& q = outi[jj];

                if (nodeId(p) == 1) {
                    spec.destruct(state(p));
                    continue;
                }

                bool allZero = true;

                for (int b = 0; b < AR; ++b) {
                    if (nodeId(p) == 0) {
                        q.branch[b] = 0;
                        continue;
                    }

                    nodeid f(i, j);
                    spec.get_copy(tmpState, state(p));
                    int kk = downTable(f, b, i - 1);
                    int ii = downSpec(tmpState, i, b, kk);

                    while (ii != 0 && kk != 0 && ii != kk) {
                        if (ii < kk) {
                            assert(kk >= 1);
                            kk = downTable(f, 0, ii);
                        }
                        else {
                            assert(ii >= 1);
                            ii = downSpec(tmpState, ii, 0, kk);
                        }
                    }

                    if (ii <= 0 || kk <= 0) {
                        if (ii == 0 || kk == 0) {
                            q.branch[b] = 0;
                        }
                        else {
                            if (oneSrcPtr.empty()) { // the first 1-terminal candidate
                                spec.get_copy(one, tmpState);
                                q.branch[b] = 1;
                                oneSrcPtr.push_back(tdzdd::NodeBranchId(i, jj, b));
                            }
                            else {
                                switch (spec.merge_states(one, tmpState)) {
                                case 1:
                                    while (!oneSrcPtr.empty()) {
                                        tdzdd::NodeBranchId const& nbi =
                                                oneSrcPtr.back();
                                        assert(nbi.row >= i);
                                        output[nbi.row][nbi.col].branch[nbi.val] =
                                                0;
                                        oneSrcPtr.pop_back();
                                    }
                                    spec.destruct(one);
                                    spec.get_copy(one, tmpState);
                                    q.branch[b] = 1;
                                    oneSrcPtr.push_back(tdzdd::NodeBranchId(i, jj, b));
                                    break;
                                case 2:
                                    q.branch[b] = 0;
                                    break;
                                default:
                                    q.branch[b] = 1;
                                    oneSrcPtr.push_back(tdzdd::NodeBranchId(i, jj, b));
                                    break;
                                }
                            }
                            allZero = false;
                        }
                    }
                    else {
                        assert(ii == f.row() && ii == kk && ii < i);
                        if (work[ii].empty()) work[ii].resize(input[ii].size());
                        SpecNode* pp = work[ii][f.col()].alloc_front(pools[ii],
                                specNodeSize);
                        spec.get_copy(state(pp), tmpState);
                        srcPtr(pp) = &q.branch[b];
                        if (ii < lowestChild) lowestChild = ii;
                        allZero = false;
                    }

                    spec.destruct(tmpState);
                }

                spec.destruct(state(p));
                ++jj;
                if (allZero) ++deadCount;
            }
        }

        work[i].clear();
        pools[i].clear();
        spec.destructLevel(i);
        sweeper.update(i, lowestChild, deadCount);
    }

private:
    int downTable(nodeid& f, int b, int zerosupLevel) const {
        if (zerosupLevel < 0) zerosupLevel = 0;

        f = input.child(f, b);
        while (f.row() > zerosupLevel) {
            f = input.child(f, 0);
        }
        return (f == 1) ? -1 : f.row();
    }

    int downSpec(void* p, int level, int b, int zerosupLevel) {
        if (zerosupLevel < 0) zerosupLevel = 0;
        assert(level > zerosupLevel);

        int i = spec.get_child(p, level, b);
        while (i > zerosupLevel) {
            i = spec.get_child(p, i, 0);
        }
        return i;
    }
};


/**
 * DD dumper.
 * A node table is printed in Graphviz (dot) format.
 */
template<typename T,typename S>
class DdDumper {
    typedef S Spec;
    static int const AR = Spec::ARITY;
    static int const headerSize = 1;

    /* SpecNode
     * ┌────────┬────────┬────────┬─────
     * │ nodeId │state[0]│state[1]│ ...
     * └────────┴────────┴────────┴─────
     */
    struct SpecNode {
        nodeid nodeId;
    };

    static nodeid& nodeId(SpecNode* p) {
        return p->nodeId;
    }

    static nodeid nodeId(SpecNode const* p) {
        return p->nodeId;
    }

    static void* state(SpecNode* p) {
        return p + headerSize;
    }

    static void const* state(SpecNode const* p) {
        return p + headerSize;
    }

    template<typename SPEC>
    struct Hasher {
        SPEC const& spec;
        int const level;

        Hasher(SPEC const& spec, int level) :
                spec(spec), level(level) {
        }

        size_t operator()(SpecNode const* p) const {
            return spec.hash_code(state(p), level);
        }

        size_t operator()(SpecNode const* p, SpecNode const* q) const {
            return spec.equal_to(state(p), state(q), level);
        }
    };

    typedef tdzdd::MyHashTable<SpecNode*,Hasher<Spec>,Hasher<Spec> > UniqTable;

    static int getSpecNodeSize(int n) {
        if (n < 0)
            throw std::runtime_error("storage size is not initialized!!!");
        return headerSize + (n + sizeof(SpecNode) - 1) / sizeof(SpecNode);
    }

    Spec spec;
    int const specNodeSize;
    char* oneState;
    nodeid oneId;

    tdzdd::MyVector<tdzdd::MyList<SpecNode> > snodeTable;
    tdzdd::MyVector<UniqTable> uniqTable;
    tdzdd::MyVector<Hasher<Spec> > hasher;

public:
    explicit DdDumper(Spec const& s) :
            spec(s),
            specNodeSize(getSpecNodeSize(spec.datasize())),
            oneState(0),
            oneId(1) {
    }

    ~DdDumper() {
        if (oneState) {
            spec.destruct(oneState);
            delete[] oneState;
        }
    }

    /**
     * Dumps the node table in Graphviz (dot) format.
     * @param os the output stream.
     * @param title title label.
     */
    void dump(std::ostream& os, std::string title) {
        if (oneState) {
            spec.destruct(oneState);
        }
        else {
            oneState = new char[spec.datasize()];
        }
        int n = spec.get_root(oneState);

        os << "digraph \"" << title << "\" {\n";

        if (n == 0) {
            if (!title.empty()) {
                os << "  labelloc=\"t\";\n";
                os << "  label=\"" << title << "\";\n";
            }
        }
        else if (n < 0) {
            os << "  \"^\" [shape=none,label=\"" << title << "\"];\n";
            os << "  \"^\" -> \"" << oneId << "\" [style=dashed" << "];\n";
            os << "  \"" << oneId << "\" ";
            os << "[shape=square,label=\"⊤\"];\n";
        }
        else {
            nodeid root(n, 0);

            for (int i = n; i >= 1; --i) {
                os << "  " << i << " [shape=none,label=\"";
                spec.printLevel(os, i);
                os << "\"];\n";
            }
            for (int i = n - 1; i >= 1; --i) {
                os << "  " << (i + 1) << " -> " << i << " [style=invis];\n";
            }

            os << "  \"^\" [shape=none,label=\"" << title << "\"];\n";
            os << "  \"^\" -> \"" << root << "\" [style=dashed" << "];\n";

            snodeTable.init(n + 1);
            SpecNode* p = snodeTable[n].alloc_front(specNodeSize);
            spec.destruct(oneState);
            spec.get_copy(state(p), oneState);
            nodeId(p) = root;

            uniqTable.clear();
            uniqTable.reserve(n + 1);
            hasher.clear();
            hasher.reserve(n + 1);
            for (int i = 0; i <= n; ++i) {
                hasher.push_back(Hasher<Spec>(spec, i));
                uniqTable.push_back(UniqTable(hasher.back(), hasher.back()));
            }

            for (int i = n; i >= 1; --i) {
                dumpStep(os, i);
            }

            for (size_t j = 2; j < oneId.code(); ++j) {
                os << "  \"" << nodeid(j) << "\" ";
                os << "[style=invis];\n";
            }
            os << "  \"" << oneId << "\" ";
            os << "[shape=square,label=\"⊤\"];\n";
        }

        os << "}\n";
        os.flush();
    }

private:
    void dumpStep(std::ostream& os, int i) {
        tdzdd::MyList<SpecNode> &snodes = snodeTable[i];
        size_t const m = snodes.size();
        tdzdd::MyVector<char> tmp(spec.datasize());
        void* const tmpState = tmp.data();
        tdzdd::MyVector<Node<T> > nodeList(m);

        for (size_t j = m - 1; j + 1 > 0; --j, snodes.pop_front()) {
            nodeid f(i, j);
            assert(!snodes.empty());
            SpecNode* p = snodes.front();

            os << "  \"" << f << "\" [label=\"";
            spec.print_state(os, state(p), i);
            os << "\"];\n";

            for (int b = 0; b < AR; ++b) {
                nodeid& child = nodeList[j].branch[b];

                if (nodeId(p) == 0) {
                    child = 0;
                    continue;
                }

                spec.get_copy(tmpState, state(p));
                int ii = spec.get_child(tmpState, i, b);

                if (ii == 0) {
                    child = 0;
                }
                else if (ii < 0) {
                    if (oneId == 1) { // the first 1-terminal candidate
                        oneId = 2;
                        spec.destruct(oneState);
                        spec.get_copy(oneState, tmpState);
                        child = oneId;
                    }
                    else {
                        switch (spec.merge_states(oneState, tmpState)) {
                        case 1:
                            oneId = oneId.code() + 1;
                            spec.destruct(oneState);
                            spec.get_copy(oneState, tmpState);
                            child = oneId;
                            break;
                        case 2:
                            child = 0;
                            break;
                        default:
                            child = oneId;
                            break;
                        }
                    }
                }
                else {
                    SpecNode* pp = snodeTable[ii].alloc_front(specNodeSize);
                    size_t jj = snodeTable[ii].size() - 1;
                    spec.get_copy(state(pp), tmpState);

                    SpecNode*& pp0 = uniqTable[ii].add(pp);
                    if (pp0 == pp) {
                        nodeId(pp) = child = nodeid(ii, jj);
                    }
                    else {
                        switch (spec.merge_states(state(pp0), state(pp))) {
                        case 1:
                            nodeId(pp0) = 0;
                            nodeId(pp) = child = nodeid(ii, jj);
                            pp0 = pp;
                            break;
                        case 2:
                            child = 0;
                            spec.destruct(state(pp));
                            snodeTable[ii].pop_front();
                            break;
                        default:
                            child = nodeId(pp0);
                            spec.destruct(state(pp));
                            snodeTable[ii].pop_front();
                            break;
                        }
                    }
                }

                spec.destruct(tmpState);
            }

            spec.destruct(state(p));
        }

        for (size_t j = 0; j < m; ++j) {
            for (int b = 0; b < AR; ++b) {
                nodeid f(i, j);
                nodeid child = nodeList[j].branch[b];
                if (child == 0) continue;
                if (child == 1) child = oneId;

                os << "  \"" << f << "\" -> \"" << child << "\"";

                os << " [style=";
                if (b == 0) {
                    os << "dashed";
                }
                else {
                    os << "solid";
                    if (AR > 2) {
                        os << ",color="
                           << ((b == 1) ? "blue" : (b == 2) ? "red" : "green");
                    }
                }
                os << "];\n";
            }
        }

        os << "  {rank=same; " << i;
        for (size_t j = 0; j < m; ++j) {
            os << "; \"" << nodeid(i, j) << "\"";
        }
        os << "}\n";

        uniqTable[i - 1].clear();
        spec.destructLevel(i);
    }
};
