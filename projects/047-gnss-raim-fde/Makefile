CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic

all: run_raim test_raim

run_raim: run_raim.cpp raim.cpp raim.hpp
	$(CXX) $(CXXFLAGS) -o $@ run_raim.cpp raim.cpp

test_raim: test_raim.cpp raim.cpp raim.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_raim.cpp raim.cpp

test: test_raim
	./test_raim
	python3 verify_stats.py ./test_raim

run: run_raim
	./run_raim

clean:
	rm -f run_raim test_raim

.PHONY: all test run clean
