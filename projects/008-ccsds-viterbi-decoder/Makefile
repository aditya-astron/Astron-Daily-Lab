CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic

SRC = viterbi.c channel.c
HDR = viterbi.h channel.h

all: ber_sim test_viterbi

ber_sim: ber_sim.c $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ ber_sim.c $(SRC) -lm

test_viterbi: test_viterbi.c $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ test_viterbi.c $(SRC) -lm

test: test_viterbi
	./test_viterbi

run: ber_sim
	./ber_sim

clean:
	rm -f ber_sim test_viterbi

.PHONY: all test run clean
