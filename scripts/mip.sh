#!/bin/bash

set -e
#get the root dir (1st ancestor of the location where this script is stored)
SRC_DIR=`dirname "$BASH_SOURCE"`/..

# Invoke MIP on the rdlp instance files in paralel. All output is collected in a single output file
# Only output lines marked with a 'tag' are preserved in the final output. This is convenient to remove debug output.
function runBenchmark(){
    instances=($(ls /home/daniel/Dropbox/PapersOR/PM/implementation/instances/wt{40,50,100,150}/*{1,6}.dat))
    machines=(2 4)
    solvers=($(seq 10))
    printf "%s\n" "${instances[@]}" | xargs -I{} printf "{} %s\n" "${machines[@]}" | xargs -I{} printf "-a %s {}\n" "${solvers[@]}" | parallel --no-notice -P 6 --eta --colsep ' ' "build/PMRelease -f 4 -S 1 {}";
}


#switch to the root directory. This allows us to invoke this script from any directory. Then run the benchmark.
pushd $SRC_DIR
runBenchmark
popd
