#!/bin/bash

echo "NOTE: this script relies on the binaries blur and blur_par to exist"
make clean > /dev/null 2>&1
make > /dev/null 2>&1

status=0
red=$(tput setaf 1)
reset=$(tput sgr0)

./blur_opt 15 "data/im1.ppm" "./data_o/blur_im1_par.ppm"
./blur_opt 15 "data/im2.ppm" "./data_o/blur_im2_par.ppm"
./blur_opt 15 "data/im3.ppm" "./data_o/blur_im3_par.ppm"
./blur_opt 15 "data/im4.ppm" "./data_o/blur_im4_par.ppm"

./blur 15 "data/im1.ppm" "./data_o/blur_im1.ppm"
./blur 15 "data/im2.ppm" "./data_o/blur_im2.ppm"
./blur 15 "data/im3.ppm" "./data_o/blur_im3.ppm"
./blur 15 "data/im4.ppm" "./data_o/blur_im4.ppm"

for image in im1 im2 im3 im4
do
    
    if ! cmp -s "./data_o/blur_${image}.ppm" "./data_o/blur_${image}_par.ppm"
    then
        echo "${red}Error: Incongruent output data detected when blurring image $image.ppm with ${reset}"
        status=1
    fi

    rm "./data_o/blur_${image}_par.ppm"
done


exit $status