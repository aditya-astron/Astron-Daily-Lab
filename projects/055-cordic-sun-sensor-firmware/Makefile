# Astron Built this AP
# Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wconversion -Wno-sign-conversion -D_DEFAULT_SOURCE
FW = cordic.c sun_sensor.c

all: test_cordic run_study

test_cordic: test_cordic.c $(FW) cordic.h sun_sensor.h common.h
	$(CC) $(CFLAGS) -o $@ test_cordic.c $(FW) -lm

run_study: run_study.c $(FW) cordic.h sun_sensor.h common.h
	$(CC) $(CFLAGS) -o $@ run_study.c $(FW) -lm

test: test_cordic
	./test_cordic

run: run_study
	./run_study

clean:
	rm -f test_cordic run_study

.PHONY: all test run clean
