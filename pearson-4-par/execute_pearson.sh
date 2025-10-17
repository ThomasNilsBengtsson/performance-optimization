#!/bin/bash

echo "Making..."
make clean > /dev/null 2>&1
make > /dev/null 2>&1

images=("128.data" "256.data" "512.data" "1024.data")

for thread in 16
do
    for img in "${images[@]}"; do
        output="output_${img}"
        echo "Running pearson on $img..."
        for i in {1..5}; do
            echo "Iteration $i and thread $thread:"
            /usr/bin/time -v ./pearson_par "data/$img" "it4_simd/pearsonOutput/$output" $thread 2>&1 #| grep -E "Percent of CPU this job got|Elapsed \(wall clock\)|Maximum resident set size|Major \(requiring I/O\) page faults|Minor \(reclaiming a frame\) page faults|Voluntary context switches|Involuntary context switches|File system inputs|File system outputs|Page size|Exit status"
        done
        echo "-----------------------------------------"
    done
done
