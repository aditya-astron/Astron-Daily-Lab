CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11 -D_DEFAULT_SOURCE
LDLIBS = -lm
OBJ = pi_q16.o wheel.o loop.o

all: wheel_demo test_loop

wheel_demo: wheel_demo.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test_loop: test_loop.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test: test_loop
	./test_loop

clean:
	rm -f *.o wheel_demo test_loop

.PHONY: all test clean
