#include <docopt/docopt.h>
#include <fmt/core.h>
#include "Parms.h"
#include "wctprivate.h"

int main(int ac, const char** av) {
    int val = 0;
    try {
        Problem problem(ac, av);
    } catch (docopt::DocoptExitHelp const&) {
        fmt::print(USAGE);
    } catch (docopt::DocoptExitVersion const&) {
        fmt::print("PM 0.1\n");
    } catch (docopt::DocoptLanguageError const& error) {
        fmt::print(stderr, "Docopt usage string could not be parsed\n");
        fmt::print(stderr, error.what());
    } catch (docopt::DocoptArgumentError const& error) {
        fmt::print("{}\n", error.what());
        fmt::print(USAGE);
    } catch (std::exception& e) {
        fmt::print(stderr, "{}", e.what());
    } catch (...) {
        fmt::print(stderr, "error: unknown exceptions");
    }
    return val;
}
