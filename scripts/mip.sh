#!/bin/bash

set -e
#get the root dir (1st ancestor of the location where this script is stored)
SRC_DIR=`dirname "$BASH_SOURCE"`/..

# Invoke MIP on the rdlp instance files in paralel. All output is collected in a single output file
# Only output lines marked with a 'tag' are preserved in the final output. This is convenient to remove debug output.
function runBenchmark(){
    
    instances=($(ls /home/daniel/Dropbox/PapersOR/PM/implementation/instances/wt{40,50,100}/*{1,6}.dat))
    # outputdir=./output
    # filename=MIPBenchmark.txt
    

    # message="MIP benchmark, 5 min, 4 threads"
    # header="name	nrCustomers    nrLocations LB   UB	feasible    optimal runtime"
    # prefix="tag: "
    
    #test whether output dir exists
    # if ! [ -d "$outputdir" ]; then
    #     mkdir ${outputdir}
    # fi

    #create the output file and write the message and header to this file
    # outputfile=$outputdir/$filename
    # touch $outputfile
    # echo "${prefix}${message}" > $outputfile
    # echo "${prefix}${header}" >> $outputfile
    
    #run the benchmark by invoking instances in parallel. -P specifies the number of instances solved in parallel.
    #additional parameters can be specified: printf "param1 param2 %s\n" "${instances[@]}"
    printf "%s 4\n" "${instances[@]}" | parallel --no-notice -P 4 --eta --colsep ' ' "build/PMRelease -a 4 -S 2 {}"
    

    #only preserve tagged lines, and remove the tag
    # mv $outputfile "${outputfile}_tmp"
    # grep ${prefix} "${outputfile}_tmp" | sed -e "s/$prefix//g" > $outputfile
    # rm "${outputfile}_tmp"

}


#switch to the root directory. This allows us to invoke this script from any directory. Then run the benchmark.
pushd $SRC_DIR
runBenchmark
popd
