#include "ZeroHalfCuts.hpp"
#include <memory>
#include <vector>
#include "ModelInterface.hpp"
#include "NodeBddTable.hpp"
#include "NodeId.hpp"

ZeroHalfCuts::ZeroHalfCuts(std::vector<BddCoeff>* _lp_sol,
                           ReformulationModel* _rmp_model, NodeId _root,
                           NodeTableEntity<>* _table)
    : env(std::make_unique<GRBEnv>()),
      model(std::make_unique<GRBModel>(*env)),
      lp_sol(_lp_sol),
      rmp_model(_rmp_model),
      root(_root),
      table(_table) {
    generate_model();
}

bool ZeroHalfCuts::add_cuts() {}
void ZeroHalfCuts::generate_model() {}

void ZeroHalfCuts::generate_cuts() {}

void ZeroHalfCuts::init_dfs() {}
