# -*- makefile -*-

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
EXT = linux
GCC_FLAGS1 = -fPIC -Wl,-Bsymbolic-functions -c -O3
GCC_FLAGS2 = -shared -O3 -Wl,-Bsymbolic-functions,-soname,transit_linux.so
endif
ifeq ($(UNAME_S),Darwin)
EXT = mac
GCC_FLAGS1 = -fPIC -c
GCC_FLAGS2 = -shared -Wl,-install_name,transit_mac.so
endif

GCC = gcc

.PHONY: all
.SILENT: all

all:
	echo "Compiling C source code..."
	${GCC} ${GCC_FLAGS1} pysyzygy/transit.c
	echo "Generating shared library..."
	gcc ${GCC_FLAGS2} -o transit_${EXT}.so transit.o -lc
	rm transit.o
	mv transit_${EXT}.so pysyzygy/.
	echo "Install successful."