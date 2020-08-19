file(GLOB SOURCE_FILES
          "src/*.cc"
          "src/*.cpp"
          src/alloc.c
          src/binomial-heap.c
          src/branch_and_bound.c
          src/conflict_branching.c
          src/copy.c src/graph.c
          src/greedy.c
          src/heap.c
          src/interval.c
          src/io.c
          src/job.c
          src/ksubset.c
          src/localsearch.c
          src/lowerbound.c
          src/lp.c
          src/model.c
          src/partlist.c
          src/preprocess.c
          src/scheduleset.c
          src/solution.c
          src/sortus.c
          src/time_bis.c
          src/time.c
          src/util.c
          src/wct.c
          src/wctparms.c
)

set(sources
    ${SOURCE_FILES}
)

set(exe_sources
		src/mainwct.c
		${sources}
)

set(headers
    ${PROJECT_SOURCE_DIR}/include
)
