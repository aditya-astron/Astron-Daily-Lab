CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

all: geo_longitude verify

geo_longitude: src/geo_longitude.cpp include/gravity.hpp include/geo.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

verify: tests/verify.cpp include/gravity.hpp include/geo.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

run: geo_longitude
	./geo_longitude

test: verify
	./verify

clean:
	rm -f geo_longitude verify

.PHONY: all run test clean
