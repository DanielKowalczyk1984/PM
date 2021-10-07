#ifndef __CARDINALITYPATHS_H__
#define __CARDINALITYPATHS_H__
#include <array>                             // for array
#include <boost/multiprecision/cpp_int.hpp>  // for cpp_int
#include <cstddef>                           // for size_t
#include "ModernDD/NodeBddEval.hpp"          // for Eval
#include "ModernDD/NodeBddTable.hpp"         // for NodeTableEntity
#include "ModernDD/NodeId.hpp"               // for NodeId
#include "NodeBdd.hpp"                       // for NodeBdd

class CardinalityPaths : public Eval<NodeBdd, boost::multiprecision::cpp_int> {
    using cpp_int = boost::multiprecision::cpp_int;

   public:
    CardinalityPaths() = default;

    cpp_int get_objective(NodeBdd& n) const override {
        return n.get_nb_paths();
    }

    void initialize_root_node(NodeBdd& n) const override { n.reset_nb_paths(); }

    void initialize_node(NodeBdd& n) const override { n.reset_nb_paths(); }

    void evalNode(NodeBdd& n) const override {
        auto table_tmp = Eval<NodeBdd, cpp_int>::get_table();
        if (n[0] == 1) {
            n.update_nb_paths();
        } else if (n[0] > 1) {
            n.update_nb_paths(table_tmp->node(n[0]).get_nb_paths());
        }

        if (n[1] == 1) {
            n.update_nb_paths();
        } else {
            n.update_nb_paths(table_tmp->node(n[1]).get_nb_paths());
        }
    }
};

class BackwardDistance : public Eval<NodeBdd, std::array<int, 2>> {
   public:
    BackwardDistance() = default;

    std::array<int, 2> get_objective(
        [[maybe_unused]] NodeBdd& n) const override {
        return {0, 0};
    }

    void initialize_root_node([[maybe_unused]] NodeBdd& n) const override {}

    void initialize_node([[maybe_unused]] NodeBdd& n) const override {}

    void evalNode(NodeBdd& n) const override {
        auto table_tmp = Eval<NodeBdd, std::array<int, 2>>::get_table();
        for (size_t i = 0UL; i < 2; ++i) {
            auto& cur_node = table_tmp->node(n[i]);
            n.update_backward_distance(cur_node.get_backward_distance(), i);
        }
    }
};
#endif  // __CARDINALITYPATHS_H__