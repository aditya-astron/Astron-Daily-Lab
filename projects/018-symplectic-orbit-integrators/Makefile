CXX ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -std=c++17

all: orbit_demo test_integrators

orbit_demo: main.o integrators.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_integrators: test_integrators.o integrators.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp integrators.hpp
	$(CXX) $(CXXFLAGS) -c $<

test: test_integrators
	./test_integrators

clean:
	rm -f *.o orbit_demo test_integrators

.PHONY: all test clean
