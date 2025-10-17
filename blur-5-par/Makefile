# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++-13
CXXFLAGS=-std=c++17 -g -Wunused -Wall -Wunused -O3 -march=native -fopenmp #-funroll-loops -fopt-info-loop-optimized #-fopt-info-vec-all

all: blur_par

blur_par: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur_par

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o  -fopt-info-vec-all

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o 

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

clean:
	rm -rf blur_par *.ppm *.o *.dSYM 2> /dev/null
