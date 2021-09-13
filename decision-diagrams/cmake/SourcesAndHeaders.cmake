set(sources
    src/tmp.cpp
)

set(exe_sources
		src/main.cpp
		${sources}
)

set(headers
    include/ModernDD/NodeBase.hpp  
    include/ModernDD/NodeBddBuilder.hpp
    include/ModernDD/NodeBddDumper.hpp
    include/ModernDD/NodeBddEval.hpp
    include/ModernDD/NodeBddReducer.hpp
    include/ModernDD/NodeBddSpec.hpp
    include/ModernDD/NodeBddStructure.hpp
    include/ModernDD/NodeBddSweeper.hpp
    include/ModernDD/NodeBddTable.hpp
    include/ModernDD/NodeBranchId.hpp
    include/ModernDD/NodeId.hpp
)

set(test_sources
  src/tmp_test.cpp
  src/example1.cpp
  src/example2.cpp
  src/test.cpp
  src/testRandomDd.cpp
  src/testSizeConstraint.cpp
)
