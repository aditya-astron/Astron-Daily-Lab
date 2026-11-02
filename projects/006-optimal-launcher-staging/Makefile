CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic

all: staging test_staging

staging: main.cpp staging.cpp staging.hpp
	$(CXX) $(CXXFLAGS) -o $@ main.cpp staging.cpp

test_staging: test_staging.cpp staging.cpp staging.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_staging.cpp staging.cpp

test: test_staging
	./test_staging

run: staging
	./staging

clean:
	rm -f staging test_staging

.PHONY: all test run clean
