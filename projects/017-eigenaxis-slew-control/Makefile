CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11 -D_DEFAULT_SOURCE
LDLIBS = -lm

all: slew_demo test_slew

slew_demo: slew_demo.o slew.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test_slew: test_slew.o slew.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

slew.o slew_demo.o test_slew.o: slew.h scenario.h

test: test_slew
	./test_slew

clean:
	rm -f *.o slew_demo test_slew *.csv

.PHONY: all test clean
