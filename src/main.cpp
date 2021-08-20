#include <docopt/docopt.h>  // for DocoptArgumentError, DocoptLanguageError
#include <fmt/core.h>       // for print
#include <cstdio>           // for stderr
#include <exception>        // for exception
#include <string>           // for string
#include "Problem.h"        // for Problem
#include "Usage.hpp"        // for USAGE

int main(int argc, const char** argv) {
    int val = 0;
    try {
        Problem problem(argc, argv);
    } catch (docopt::DocoptExitHelp const&) {
        fmt::print("{}", USAGE);
    } catch (docopt::DocoptExitVersion const&) {
        fmt::print("PM 0.1\n");
    } catch (docopt::DocoptLanguageError const& error) {
        fmt::print(stderr, "Docopt usage string could not be parsed\n");
        fmt::print(stderr, "{}\n", error.what());
    } catch (docopt::DocoptArgumentError const& error) {
        fmt::print("{}\n", error.what());
        fmt::print("{}", USAGE);
    } catch (std::exception& e) {
        fmt::print(stderr, "{}\n", e.what());
    } catch (...) {
        fmt::print(stderr, "error: unknown exceptions\n");
    }
    return val;
}
