CC     ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Wpedantic

all: rs_link test_rs librs.so

rs_link: rs_link.c rs.c rs.h
	$(CC) $(CFLAGS) -o $@ rs_link.c rs.c -lm

test_rs: test_rs.c rs.c rs.h
	$(CC) $(CFLAGS) -o $@ test_rs.c rs.c

# shared build for the independent cross-check in verify_galois.py
librs.so: rs.c rs.h
	$(CC) $(CFLAGS) -fPIC -shared -o $@ rs.c

test: test_rs
	./test_rs

verify: librs.so
	python3 verify_galois.py

run: rs_link
	./rs_link

clean:
	rm -f rs_link test_rs librs.so

.PHONY: all test verify run clean
