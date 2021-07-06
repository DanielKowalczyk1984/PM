#ifndef __BRANCHNODE_H__
#define __BRANCHNODE_H__

#include <fmt/core.h>
#include <cstddef>
#include <limits>
#include <memory>
#include "Instance.h"
#include "PricerSolverBase.hpp"
#include "branch-and-bound/btree.h"
#include "branch-and-bound/state.h"

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

    static constexpr size_t NbCandidates = 16;

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

struct BranchCandidate {
    double score{std::numeric_limits<double>::lowest()};
    bool   empty{true};
    std::array<std::unique_ptr<NodeData>, 2> data_child_nodes;

    BranchCandidate() = default;
    BranchCandidate(const BranchCandidate&) = delete;
    BranchCandidate(BranchCandidate&&) = default;
    BranchCandidate& operator=(const BranchCandidate&) = delete;
    BranchCandidate& operator=(BranchCandidate&&) = default;
    ~BranchCandidate() = default;

    BranchCandidate(double                                     _score,
                    std::array<std::unique_ptr<NodeData>, 2>&& child_nodes);
    BranchCandidate(double _score);
};
#endif  // __BRANCHNODE_H__