CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic

all: tm_link test_ccsds

tm_link: tm_link.c ccsds.c ccsds.h
	$(CC) $(CFLAGS) -o $@ tm_link.c ccsds.c -lm

test_ccsds: test_ccsds.c ccsds.c ccsds.h
	$(CC) $(CFLAGS) -o $@ test_ccsds.c ccsds.c

test: test_ccsds
	./test_ccsds

run: tm_link
	./tm_link

clean:
	rm -f tm_link test_ccsds

.PHONY: all test run clean
