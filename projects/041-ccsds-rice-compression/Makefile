CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -D_DEFAULT_SOURCE
LDLIBS  = -lm

all: rice_demo rice_tool test_rice

rice_demo: main.c rice.c rice.h
	$(CC) $(CFLAGS) -o $@ main.c rice.c $(LDLIBS)

rice_tool: rice_tool.c rice.c rice.h
	$(CC) $(CFLAGS) -o $@ rice_tool.c rice.c $(LDLIBS)

test_rice: test_rice.c rice.c rice.h
	$(CC) $(CFLAGS) -o $@ test_rice.c rice.c $(LDLIBS)

test: test_rice rice_tool
	./test_rice
	python3 test_libaec_interop.py

run: rice_demo
	./rice_demo

clean:
	rm -f rice_demo rice_tool test_rice

.PHONY: all test run clean
