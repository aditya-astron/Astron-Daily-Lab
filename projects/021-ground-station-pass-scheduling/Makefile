CXX ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -std=c++17

all: pass_scheduler test_scheduler

pass_scheduler: main.o scheduler.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_scheduler: test_scheduler.o scheduler.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp scheduler.hpp
	$(CXX) $(CXXFLAGS) -c $<

test: test_scheduler
	./test_scheduler

clean:
	rm -f *.o pass_scheduler test_scheduler

.PHONY: all test clean
