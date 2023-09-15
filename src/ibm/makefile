# makefile to compile individual-based model 
# on iterated prisoner's dilemma written in C++

EXE=ipd_simulation.exe
CPP=simulation.cpp main.cpp individual.cpp
HPP=simulation.hpp parameters.hpp individual.hpp
CXX=g++
CXXFLAGS=-Wall -std=c++20 -O3


$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)
