// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef __INSTANCE_H__
#define __INSTANCE_H__
#include <fmt/format.h>  // for format
#include <chrono>        // for filesystem
#include <cstddef>       // for size_t
#include <exception>     // for exception
#include <filesystem>    // for path
#include <limits>        // for numeric_limits
#include <memory>        // for shared_ptr
#include <string>        // for string
#include <utility>       // for pair
#include <vector>        // for vector
struct Interval;
struct Job;
struct Parms;

namespace fs = std::filesystem;

struct Instance {
    fs::path                                path_to_instance{};
    std::vector<std::shared_ptr<Job>>       jobs;
    std::vector<std::shared_ptr<Interval>>  intervals;
    std::vector<std::pair<Job*, Interval*>> vector_pair;

    /** Summary of jobs */
    std::string pname;
    size_t      nb_jobs{};
    size_t      nb_machines{};
    int         p_sum{};
    int         pmax{};
    int         pmin{std::numeric_limits<int>::max()};
    int         dmax{};
    int         dmin{std::numeric_limits<int>::max()};
    int         H_min{};
    int         H_max{};
    int         off{};

    Instance(Parms const& _parms);

    Instance(const Instance&) = default;
    Instance& operator=(const Instance&) = delete;
    Instance& operator=(Instance&&) = delete;
    Instance(Instance&&) = default;
    ~Instance() = default;

    class InstanceException : public std::exception {
       public:
        InstanceException(const char* const msg = nullptr) : errmsg(msg) {}

        [[nodiscard]] const char* what() const noexcept override {
            return (errmsg);
        }

       private:
        const char* errmsg;
    };

   private:
    void calculate_H_max_H_min();
    void find_division();
};

template <>
struct fmt::formatter<Instance> {
    char presentation = 'v';

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && (*it == 'v' || *it == 'n'))
            presentation = *it++;

        // Check if reached the end of the range:
        if (it != end && *it != '}')
            throw format_error("invalid format");

        return it;
    }

    template <typename FormatContext>
    auto format(const Instance& inst, FormatContext& ctx)
        -> decltype(ctx.out()) {
        if (presentation == 'n') {
            return format_to(ctx.out(), "n,m,NameInstance");
        } else {
            return format_to(ctx.out(), "{},{},{}", inst.nb_jobs,
                             inst.nb_machines, inst.pname);
        }
    }
};
#endif  // __INSTANCE_H__
