# tsunami Makefile

FC = gfortran
FCFLAGS = -O3

.PHONY: all clean

all: ch02 ch03 ch04

ch02:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch02

ch03:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch03

ch04:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch04

clean:
	$(MAKE) clean --directory=src/ch02
	$(MAKE) clean --directory=src/ch03
	$(MAKE) clean --directory=src/ch04
