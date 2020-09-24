#ifndef __ZEROHALFCUTS_H__
#define __ZEROHALFCUTS_H__

#include <memory>
#include <vector>
#include "ModelInterface.hpp"
#include "NodeBdd.hpp"
#include "NodeBddTable.hpp"
#include "NodeId.hpp"
#include "gurobi_c++.h"

class ZeroHalfCuts {
   public:
    ZeroHalfCuts(int                 _nb_jobs,
                 int                 _nb_machines,
                 ReformulationModel* _rmp_model,
                 NodeId&             _root,
                 NodeTableEntity<>*  _table);
    ZeroHalfCuts(ZeroHalfCuts&&) = default;
    ZeroHalfCuts(const ZeroHalfCuts&) = delete;
    ZeroHalfCuts& operator=(ZeroHalfCuts&&) = default;
    ZeroHalfCuts& operator=(const ZeroHalfCuts&) = delete;

    ~ZeroHalfCuts() = default;
    bool                                            add_cuts();
    void                                            generate_cuts();
    std::vector<std::shared_ptr<ConstraintGeneric>> get_cut_list();

   private:
    std::unique_ptr<GRBEnv>                         env;
    std::unique_ptr<GRBModel>                       model;
    int                                             nb_jobs;
    int                                             nb_machines;
    ReformulationModel*                             rmp_model;
    NodeId                                          root;
    NodeTableEntity<>*                              table;
    std::vector<NodeId>                             node_ids{};
    std::vector<NodeId>                             node_ids_lift{};
    std::vector<GRBVar>                             jobs_var;
    GRBVar                                          q;
    std::vector<std::shared_ptr<ConstraintGeneric>> cut_list;

    void generate_model();
    void init_table();
    void init_coeff_cut();
    void init_coeff_node(NodeBdd<>& node);
    void construct_cut();
    void lift_operator();
    void dfs(const NodeId& v);
    void dfs_lift(const NodeId& v);
    void dfs_evaluate(const NodeId& v);
};

#endif  // __ZEROHALFCUTS_H__
