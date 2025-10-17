#!/bin/bash

echo "Making..."
make clean > /dev/null 2>&1
make > /dev/null 2>&1

images=("im1.ppm" "im2.ppm" "im3.ppm" "im4.ppm")

for img in "${images[@]}"; do
    output="output_${img}"
    echo "Running blur on $img..."
    echo "Iteration $i:"
    valgrind --tool=callgrind ./blur_opt 15 "data/$img" "./it5-simd/kcachegrind/$output" 2>&1 #| grep -E "Percent of CPU this job got|Elapsed \(wall clock\)|Maximum resident set size|Major \(requiring I/O\) page faults|Minor \(reclaiming a frame\) page faults|Voluntary context switches|Involuntary context switches|File system inputs|File system outputs|Page size|Exit status"
    echo "-----------------------------------------"
done
