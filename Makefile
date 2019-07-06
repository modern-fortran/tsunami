# tsunami Makefile

FC = gfortran
FCFLAGS = -O3

.PHONY: all clean

all: ch02 ch03

ch02:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch02

ch03:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch03

clean:
	$(MAKE) clean --directory=src/ch02
	$(MAKE) clean --directory=src/ch03
