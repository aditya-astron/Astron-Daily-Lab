CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

all: star_id test_star_id

star_id: main.cpp star_id.hpp
	$(CXX) $(CXXFLAGS) -o $@ main.cpp

test_star_id: test_star_id.cpp star_id.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_star_id.cpp

test: test_star_id
	./test_star_id

run: star_id
	./star_id

clean:
	rm -f star_id test_star_id

.PHONY: all test run clean
