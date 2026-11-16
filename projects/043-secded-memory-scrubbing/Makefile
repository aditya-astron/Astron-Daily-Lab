CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -D_DEFAULT_SOURCE
LDLIBS  = -lm

all: scrub_design test_secded

scrub_design: main.c secded.c secded.h
	$(CC) $(CFLAGS) -o $@ main.c secded.c $(LDLIBS)

test_secded: test_secded.c secded.c secded.h scrub_sim.c scrub_sim.h
	$(CC) $(CFLAGS) -o $@ test_secded.c secded.c scrub_sim.c $(LDLIBS)

test: test_secded
	./test_secded

run: scrub_design
	./scrub_design

clean:
	rm -f scrub_design test_secded

.PHONY: all test run clean
