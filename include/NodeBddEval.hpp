#ifndef __NODEBDDEVAL_H__
#define __NODEBDDEVAL_H__

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

#include "NodeBddTable.hpp"  // for NodeTableEntity

/**
 * @brief Base class of the evaluators
 * Every derived class needs an implementation of the following functions:
 * - void initializenode(T&)
 * - void initializerootnode(T&)
 * - void evalNode(T&)
 * - R get_objective(T&)
 *
 * @tparam T: type of the Node
 * @tparam R: type of the Solution
 */
template <typename T, typename R = T>
class Eval {
    NodeTableEntity<T>* table;

   public:
    [[nodiscard]] virtual bool showMessages() const { return false; }

    /**
     * @brief Initialization of the node before evaluation.
     *
     * @param n Node at which we initialize
     */
    virtual void initializenode(T& n) const = 0;

    /**
     * @brief Initialization of the root node
     *
     * @param n Root node of the decision diagram
     */
    virtual void initializerootnode(T& n) const = 0;

    /**
     * @brief Evaluate the node
     *
     * @param n Node at which we evaluate
     */
    virtual void evalNode(T& n) const = 0;

    /**
     * @brief Get the objective object by backtracking
     *
     * @param n Depends on the direction of the evaluation. If backward
     * evaluation then n is the root node otherwise we start at terminal node 1.
     * @return R Solution we want to obtain.
     */
    virtual R get_objective(T& n) const = 0;

    void set_table(NodeTableEntity<T>* _table) { table = _table; }
    NodeTableEntity<T>* get_table() const { return table; }

    /** Default base constructors */
    Eval() : table(nullptr){};
    Eval(const Eval<T, R>&) = default;
    Eval(Eval<T, R>&&) noexcept = default;
    Eval<T, R>& operator=(const Eval<T, R>&) = default;
    Eval<T, R>& operator=(Eval<T, R>&&) noexcept = default;
    virtual ~Eval() = default;
};
#endif  // __NODEBDDEVAL_H__