#ifndef __CARDINALITYPATHS_H__
#define __CARDINALITYPATHS_H__
#include <boost/multiprecision/cpp_int.hpp>
#include "NodeBdd.hpp"
#include "NodeBddEval.hpp"

class CardinalityPaths
    : public Eval<NodeBdd<>, boost::multiprecision::cpp_int> {
    using cpp_int = boost::multiprecision::cpp_int;

   public:
    CardinalityPaths() = default;

    cpp_int get_objective(NodeBdd<>& n) const override { return n.nb_paths; }

    void initializerootnode(NodeBdd<>& n) const override { n.nb_paths = 0UL; }

    void initializenode(NodeBdd<>& n) const override { n.nb_paths = 0UL; }

    void evalNode(NodeBdd<>& n) const override {
        auto table_tmp = Eval<NodeBdd<>, cpp_int>::get_table();
        if (n[0] == 1) {
            n.nb_paths += 1;
        } else if (n[0] > 1) {
            n.nb_paths = table_tmp->node(n[0]).nb_paths;
        }

        if (n[1] == 1) {
            n.nb_paths += 1;
        } else {
            n.nb_paths += table_tmp->node(n[1]).nb_paths;
        }
    }
};
#endif  // __CARDINALITYPATHS_H__