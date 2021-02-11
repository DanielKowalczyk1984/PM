file(
  GLOB
  SOURCE_FILES
  "src/*.cpp"
  src/Statistics.c
  src/greedy.c
  src/interval.c
  src/io.c
  src/job.c
  src/localsearch.c
  src/lowerbound.c
  src/lp.c
  src/model.c
  src/partlist.c
  src/preprocess.c
  src/scheduleset.c
  src/solution.c
  src/wct.c
  src/wctparms.c)

set(sources ${SOURCE_FILES})

set(exe_sources src/main.cpp)

set(headers ${PROJECT_SOURCE_DIR}/include)
