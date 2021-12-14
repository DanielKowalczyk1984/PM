#!/bin/bash

set -e
#get the root dir (1st ancestor of the location where this script is stored)
SRC_DIR=$(dirname "$BASH_SOURCE")/..

# Invoke MIP on the rdlp instance files in paralel. All output is collected in a single output file
# Only output lines marked with a 'tag' are preserved in the final output. This is convenient to remove debug output.
function runBenchmark() {
    instances=($(fd --regex "wt(040|050|100).*[16]\.dat" .))
    machines=(4 2)
    solvers=4
    printf "%s\n" "${machines[@]}" | xargs -I{} printf "%s {}\n" "${instances[@]}" | parallel --no-notice -P 25 --eta --colsep ' ' "build/PM --json settings/default_settings.json {}"
}

#switch to the root directory. This allows us to invoke this script from any directory. Then run the benchmark.
pushd $SRC_DIR
runBenchmark
popd
