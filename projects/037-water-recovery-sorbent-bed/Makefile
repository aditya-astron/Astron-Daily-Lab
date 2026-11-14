CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic

all: bed_design test_sorbent_bed

bed_design: bed_design.cpp sorbent_bed.hpp
	$(CXX) $(CXXFLAGS) -o $@ bed_design.cpp

test_sorbent_bed: test_sorbent_bed.cpp sorbent_bed.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_sorbent_bed.cpp

run: bed_design
	./bed_design

test: test_sorbent_bed
	./test_sorbent_bed

clean:
	rm -f bed_design test_sorbent_bed

.PHONY: all run test clean
