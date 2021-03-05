set(SOURCE_FILES
    src/BranchBoundTree.cpp
    src/BranchNode.cpp
    src/greedy.cpp
    src/Instance.cpp
    src/interval.cpp
    src/Interval_new.cpp
    src/io.cpp
    src/job.c
    src/localsearch.c
    src/LocalSearch_new.cpp
    src/lowerbound.cpp
    src/model.cpp
    src/ModelInterface.cpp
    src/parms.cpp
    src/partlist.c
    src/preprocess.cpp
    src/PricerSolverArcTimeDP.cpp
    src/PricerSolverBase.cpp
    src/PricerSolverBdd.cpp
    src/PricerSolverBddBackward.cpp
    src/PricerSolverBddForward.cpp
    src/PricerSolverSimpleDP.cpp
    src/PricerSolverWrappers.cpp
    src/PricerSolverZdd.cpp
    src/PricerSolverZddBackward.cpp
    src/PricerSolverZddForward.cpp
    src/PricingStabilization.cpp
    src/scheduleset.cpp
    src/SeperationSolver.cpp
    src/solution.c
    src/Solution_new.cpp
    src/StabilizationWrappers.cpp
    src/Statistics.cpp
    src/wct.cpp
    src/wctprivate.cpp
    src/ZeroHalfCuts.cpp)

set(sources ${SOURCE_FILES})

set(exe_sources src/main.cpp)

set(headers ${PROJECT_SOURCE_DIR}/include)
