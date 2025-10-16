# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++-13
CXXFLAGS=-std=c++17 -g -Wunused -Wall -Wunused -Ofast -mavx2

all: pearson_opt verify

pearson_opt: vector dataset analysis pearson.cpp
	$(CXX) $(CXXFLAGS) pearson.cpp vector.o dataset.o analysis.o -o pearson_opt

analysis: vector analysis.hpp analysis.cpp
	$(CXX) $(CXXFLAGS) -c analysis.cpp -o analysis.o 

dataset: vector dataset.hpp dataset.cpp
	$(CXX) $(CXXFLAGS) -c dataset.cpp -o dataset.o

vector: vector.hpp vector.cpp
	$(CXX) $(CXXFLAGS) -c vector.cpp -o vector.o -fopt-info-vec-all

verify: verify.c
	$(CC) verify.c -o verify

clean:
	rm -rf verify pearson_opt *.o *.dSYM 2> /dev/null
