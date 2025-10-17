/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>
#include <omp.h>

int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    int num_threads = std::atoi(argv[3]);
    omp_set_num_threads(num_threads);

    auto datasets { Dataset::read(argv[1]) };
    auto corrs { Analysis::correlation_coefficients(datasets) };
    Dataset::write(corrs, argv[2]);

    return 0;
}