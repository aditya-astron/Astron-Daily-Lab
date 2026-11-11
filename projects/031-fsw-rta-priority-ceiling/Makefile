CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic
LDLIBS = -lm

all: rta_report test_sched

rta_report: main.c sched.c fsw_tasks.c sched.h fsw_tasks.h
	$(CC) $(CFLAGS) -o $@ main.c sched.c fsw_tasks.c $(LDLIBS)

test_sched: test_sched.c sched.c fsw_tasks.c sched.h fsw_tasks.h
	$(CC) $(CFLAGS) -o $@ test_sched.c sched.c fsw_tasks.c $(LDLIBS)

test: test_sched
	./test_sched

run: rta_report
	./rta_report

clean:
	rm -f rta_report test_sched

.PHONY: all test run clean
