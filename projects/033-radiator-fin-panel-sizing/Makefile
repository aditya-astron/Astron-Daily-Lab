CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic

all: radiator_design test_radiator

radiator_design: radiator_design.cpp radiator.hpp
	$(CXX) $(CXXFLAGS) -o $@ radiator_design.cpp

test_radiator: test_radiator.cpp radiator.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_radiator.cpp

run: radiator_design
	./radiator_design

test: test_radiator
	./test_radiator

clean:
	rm -f radiator_design test_radiator

.PHONY: all run test clean
