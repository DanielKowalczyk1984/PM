// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//
// Created by daniel on 10/6/21.
//

#ifndef PM_VARIABLEKEYBASE_HPP
#define PM_VARIABLEKEYBASE_HPP

#include <boost/container_hash/extensions.hpp>  // for hash_combine
#include <cstddef>
class VariableKeyBase {
   private:
    size_t j{};
    size_t t{};
    bool   high{false};
    bool   root{false};

   public:
    VariableKeyBase(size_t _j,
                    size_t _t,
                    bool   _high = true,
                    bool   _root = false);

    VariableKeyBase();

    [[nodiscard]] size_t get_j() const;
    [[nodiscard]] size_t get_t() const;
    [[nodiscard]] bool   get_root() const;
    [[nodiscard]] bool   get_high() const;

    void set_j(size_t _j);
    void set_t(size_t _t);
    void set_root(bool _root);
    void set_high(bool _high);

    VariableKeyBase(const VariableKeyBase&);
    VariableKeyBase(VariableKeyBase&&) noexcept;
    VariableKeyBase& operator=(const VariableKeyBase&);
    VariableKeyBase& operator=(VariableKeyBase&&) noexcept;
    virtual ~VariableKeyBase() = default;

    bool operator==(const VariableKeyBase& other) const;
};

namespace std {
template <>
struct hash<VariableKeyBase> {
    std::size_t operator()(auto const& s) const noexcept {
        std::size_t seed = 0;
        boost::hash_combine(seed, s.get_j());
        boost::hash_combine(seed, s.get_t());
        boost::hash_combine(seed, s.get_high());

        return seed;  // or use boost::hash_combine
    }
};
}  // namespace std
#endif  // PM_VARIABLEKEYBASE_HPP
