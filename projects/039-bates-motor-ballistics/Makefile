CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -D_DEFAULT_SOURCE
LDLIBS  = -lm

all: motor test_ballistics

motor: main.c ballistics.c ballistics.h motor_design.h
	$(CC) $(CFLAGS) -o $@ main.c ballistics.c $(LDLIBS)

test_ballistics: test_ballistics.c ballistics.c ballistics.h motor_design.h
	$(CC) $(CFLAGS) -o $@ test_ballistics.c ballistics.c $(LDLIBS)

test: test_ballistics
	./test_ballistics

run: motor
	./motor

clean:
	rm -f motor test_ballistics

.PHONY: all test run clean
