#ifndef __ZEROHALFCUTS_H__
#define __ZEROHALFCUTS_H__

#include <gurobi_c++.h>               // for GRBVar, GRBEnv, GRBModel
#include <cstddef>                    // for size_t
#include <memory>                     // for shared_ptr, unique_ptr
#include <vector>                     // for vector
#include "ModernDD/NodeBddTable.hpp"  // for NodeTableEntity
#include "ModernDD/NodeId.hpp"        // for NodeId
#include "NodeBdd.hpp"                // for NodeBdd
class ConstraintGeneric;
class ReformulationModel;
class ZeroHalfCuts {
   public:
    ZeroHalfCuts(size_t                    _nb_jobs,
                 size_t                    _nb_machines,
                 ReformulationModel*       _rmp_model,
                 NodeId const&             _root,
                 NodeTableEntity<NodeBdd>* _table);
    ZeroHalfCuts(ZeroHalfCuts&&) = default;
    ZeroHalfCuts(const ZeroHalfCuts&) = delete;
    ZeroHalfCuts& operator=(ZeroHalfCuts&&) = default;
    ZeroHalfCuts& operator=(const ZeroHalfCuts&) = delete;

    ~ZeroHalfCuts() = default;
    bool                                            add_cuts();
    void                                            generate_cuts();
    std::vector<std::shared_ptr<ConstraintGeneric>> get_cut_list();

    static constexpr int ALIGN = 40;

   private:
    std::unique_ptr<GRBEnv>                         env;
    std::unique_ptr<GRBModel>                       model;
    size_t                                          nb_jobs;
    size_t                                          nb_machines;
    ReformulationModel*                             rmp_model;
    NodeId                                          root;
    NodeTableEntity<NodeBdd>*                       table;
    std::vector<NodeId>                             node_ids{};
    std::vector<NodeId>                             node_ids_lift{};
    std::vector<GRBVar>                             jobs_var;
    GRBVar                                          q;
    std::vector<std::shared_ptr<ConstraintGeneric>> cut_list;

    static constexpr double EPS_CUT = 1e-6;
    static constexpr double HALF = 2.0;
    static constexpr double TIMELIMIT = 50.0;

    void generate_model();
    void init_table();
    void init_coeff_cut();
    void init_coeff_node(NodeBdd* node);
    void construct_cut();
    void lift_operator();
    void dfs(const NodeId& v);
    void dfs_lift(const NodeId& v);
    void dfs_evaluate(const NodeId& v);
};

#endif  // __ZEROHALFCUTS_H__
