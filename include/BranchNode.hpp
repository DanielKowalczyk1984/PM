#ifndef __BRANCHNODE_H__
#define __BRANCHNODE_H__

#include <fmt/core.h>
#include <limits>
#include <memory>
#include "Instance.h"
#include "PricerSolverBase.hpp"
#include "branch-and-bound/btree.h"
#include "branch-and-bound/state.h"
// #include "wctprivate.h"
struct NodeData;
struct Instance;
class BranchNodeBase : public State {
   public:
    explicit BranchNodeBase(std::unique_ptr<NodeData> pd, bool isRoot = false);
    BranchNodeBase(BranchNodeBase&&) = default;
    BranchNodeBase(const BranchNodeBase&) = delete;
    BranchNodeBase& operator=(BranchNodeBase&&) = default;
    BranchNodeBase& operator=(const BranchNodeBase&) = delete;
    ~BranchNodeBase() override = default;
    [[nodiscard]] std::unique_ptr<State> clone() const override {
        return nullptr;
    }

    void branch(BTree* bt) override;
    void compute_bounds(BTree* bt) final;
    void assess_dominance(State* otherState) final;
    bool is_terminal_state() final;
    void apply_final_pruning_tests(BTree* bt) final;
    void update_data(double upper_bound) final;
    void print(const BTree* bt) const override;

    [[nodiscard]] NodeData*         get_data_ptr() const { return pd.get(); }
    [[nodiscard]] const Instance&   get_instance_info() const;
    [[nodiscard]] PricerSolverBase* get_pricersolver() const;

   private:
    std::unique_ptr<NodeData> pd;
    static constexpr double   ERROR = 1e-12;
    static constexpr double   IntegerTolerance = 1e-3;
};

class BranchNodeRelBranching : public BranchNodeBase {
   public:
    explicit BranchNodeRelBranching(std::unique_ptr<NodeData>,
                                    bool isRoot = false);
    BranchNodeRelBranching(BranchNodeRelBranching&&) = default;
    BranchNodeRelBranching(const BranchNodeRelBranching&) = delete;
    BranchNodeRelBranching& operator=(BranchNodeRelBranching&&) = default;
    BranchNodeRelBranching& operator=(const BranchNodeRelBranching&) = delete;
    ~BranchNodeRelBranching() override = default;

    void branch(BTree* bt) override;
};

struct BranchCand {
    double score{EPS_BRANCH};
    int    job{-1};
    int    t{-1};

    // std::unique_ptr<BranchNodeBase> left{};
    // std::unique_ptr<BranchNodeBase> right{};

    BranchCand() = default;

    BranchCand(double _score, int _job, int _t);

    static constexpr double EPS_BRANCH = 1e-4;
};
#endif  // __BRANCHNODE_H__