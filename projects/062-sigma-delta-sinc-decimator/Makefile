# Astron Built this AP
# Project 062 | Signal processing / Embedded firmware | Sinc^N decimation of a second-order sigma-delta stream
# Author: Aditya (aditya-astron)
# Date: 2026-11-24 IST
CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -D_DEFAULT_SOURCE

all: test_cic libcic.so

test_cic: test_cic.c cic.c cic.h
	$(CC) $(CFLAGS) -o $@ test_cic.c cic.c -lm

libcic.so: cic.c cic.h
	$(CC) $(CFLAGS) -fPIC -shared -o $@ cic.c -lm

test: test_cic libcic.so
	./test_cic
	python3 -m pytest -q test_decimation.py

run: libcic.so
	python3 decimation_study.py

clean:
	rm -f test_cic libcic.so

.PHONY: all test run clean
