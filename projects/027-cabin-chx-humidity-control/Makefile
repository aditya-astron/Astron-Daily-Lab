CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic

all: cabin_humidity test_chx

cabin_humidity: cabin_humidity.cpp chx.hpp
	$(CXX) $(CXXFLAGS) -o $@ cabin_humidity.cpp

test_chx: test_chx.cpp chx.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_chx.cpp

run: cabin_humidity
	./cabin_humidity

test: test_chx
	./test_chx

clean:
	rm -f cabin_humidity test_chx

.PHONY: all run test clean
