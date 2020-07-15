#ifndef __ZEROHALFCUTS_H__
#define __ZEROHALFCUTS_H__

#include <memory>
#include <vector>
#include "ModelInterface.hpp"
#include "NodeBddTable.hpp"
#include "NodeId.hpp"
#include "gurobi_c++.h"

class ZeroHalfCuts {
    ZeroHalfCuts(std::vector<BddCoeff>* _lp_sol, ReformulationModel* _rmp_model,
                 NodeId _root, NodeTableEntity<>* _table);
    ZeroHalfCuts(ZeroHalfCuts&&) = delete;
    ZeroHalfCuts(const ZeroHalfCuts&) = delete;
    ZeroHalfCuts& operator=(ZeroHalfCuts&&) = delete;
    ZeroHalfCuts& operator=(const ZeroHalfCuts&) = delete;

    ~ZeroHalfCuts();
    bool add_cuts();

   private:
    std::unique_ptr<GRBEnv>   env;
    std::unique_ptr<GRBModel> model;
    std::vector<BddCoeff>*    lp_sol;
    ReformulationModel*       rmp_model;
    NodeId                    root;
    NodeTableEntity<>*        table;

    void generate_model();
    void generate_cuts();
    void init_dfs();
};

#endif  // __ZEROHALFCUTS_H__
