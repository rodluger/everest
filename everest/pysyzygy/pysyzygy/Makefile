# -*- makefile -*-

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
GCC_FLAGS1 = -fPIC -Wl,-Bsymbolic-functions -c -O3
GCC_FLAGS2 = -shared -O3 -Wl,-Bsymbolic-functions,-soname,transitlib.so
endif
ifeq ($(UNAME_S),Darwin)
GCC_FLAGS1 = -fPIC -c
GCC_FLAGS2 = -shared -Wl,-install_name,transitlib.so
endif

GCC = gcc

.PHONY: all
.SILENT: all

all:
	echo "[pysyzygy] Compiling C source code..."
	${GCC} ${GCC_FLAGS1} transit.c
	echo "[pysyzygy] Generating shared library..."
	gcc ${GCC_FLAGS2} -o transitlib.so transit.o -lc
	rm transit.o
	echo "[pysyzygy] Install successful."