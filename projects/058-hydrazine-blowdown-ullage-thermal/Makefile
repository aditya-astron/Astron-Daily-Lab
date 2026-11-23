CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic

all: run_study test_blowdown

run_study: run_study.cpp blowdown.hpp design.hpp
	$(CXX) $(CXXFLAGS) -o $@ run_study.cpp

test_blowdown: test_blowdown.cpp blowdown.hpp design.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_blowdown.cpp

test: test_blowdown
	./test_blowdown

run: run_study
	./run_study

clean:
	rm -f run_study test_blowdown
