file(GLOB SOURCE_FILES
          "src/*.cpp"
          src/branch_and_bound.c
          src/conflict_branching.c
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
          src/wctparms.c
)

set(sources
    ${SOURCE_FILES}
)

set(exe_sources
		src/mainwct.c
)

set(headers
    ${PROJECT_SOURCE_DIR}/include
)
