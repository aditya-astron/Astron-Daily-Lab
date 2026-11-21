CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic

all: run_soc test_soc

run_soc: run_soc.cpp soc.hpp leo.hpp
	$(CXX) $(CXXFLAGS) -o $@ run_soc.cpp

test_soc: test_soc.cpp soc.hpp leo.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_soc.cpp

test: test_soc
	./test_soc

run: run_soc
	./run_soc

clean:
	rm -f run_soc test_soc
