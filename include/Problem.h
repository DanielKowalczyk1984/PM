#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H
#include <array>         // for array
#include <cstddef>       // for size_t
#include <exception>     // for exception
#include <functional>    // for function
#include <memory>        // for unique_ptr, shared_ptr
#include <string>        // for string
#include <vector>        // for vector
#include "Instance.h"    // for Instance
#include "Parms.h"       // for Parms
#include "Solution.hpp"  // for Sol
#include "Statistics.h"  // for Statistics
#include "lp.h"          // for wctlp
class BranchBoundTree;   // lines 21-21
class PricingStabilizationBase;
struct NodeData;  // lines 19-19
struct PricerSolverBase;
struct Column;  // lines 18-18

/**
 * wct data types nodes of branch and bound tree
 */
/**
 *  CONSTANTS NODEDATA STRUCTURE
 *
 */

enum problem_status {
    no_sol = 0,
    lp_feasible = 1,
    feasible = 2,
    meta_heuristic = 3,
    optimal = 4
};

/**
 * problem data
 */
class Problem {
   private:
    /** Different Parameters */
    Parms parms;
    /*Cpu time measurement + Statistics*/
    Statistics stat;
    /** Instance data*/
    Instance instance;

    std::unique_ptr<BranchBoundTree> tree;
    std::unique_ptr<NodeData>        root_pd;

    problem_status status;

    /* Best Solution*/
    Sol opt_sol;

    static constexpr auto EPS = 1e-6;

   public:
    /** All methods of problem class */
    int  to_screen();
    void to_csv();
    void solve();
    /** Heuristic related */
    void heuristic();
    /** Constructors */
    Problem(int argc, const char** argv);
    Problem(const Problem&) = delete;
    Problem(Problem&&) = delete;
    Problem& operator=(const Problem&) = delete;
    Problem& operator=(Problem&&) = delete;
    ~Problem();
    friend NodeData;

    class ProblemException : public std::exception {
       public:
        ProblemException(const char* const msg = nullptr) : errmsg(msg) {}

        [[nodiscard]] const char* what() const noexcept override {
            return (errmsg);
        }

       private:
        const char* errmsg;
    };
};

#endif
