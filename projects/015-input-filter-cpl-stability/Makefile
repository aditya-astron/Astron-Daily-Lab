CXX ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -std=c++17

all: filter_demo test_filter

filter_demo: main.o filter.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_filter: test_filter.o filter.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp filter.hpp
	$(CXX) $(CXXFLAGS) -c $<

test: test_filter
	./test_filter

clean:
	rm -f *.o filter_demo test_filter

.PHONY: all test clean
