# tsunami Makefile

FC = gfortran
FCFLAGS = -O3

.PHONY: all clean

all: ch02

ch02:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch02

clean:
	$(MAKE) clean --directory=src/ch02
