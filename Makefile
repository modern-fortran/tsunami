# tsunami Makefile

FC = caf
FCFLAGS = -O3

.PHONY: all clean

all: ch02 ch03 ch04 ch07 ch10

ch02:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch02

ch03:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch03

ch04:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch04

ch07:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch07

ch10:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/ch10

clean:
	$(MAKE) clean --directory=src/ch02
	$(MAKE) clean --directory=src/ch03
	$(MAKE) clean --directory=src/ch04
	$(MAKE) clean --directory=src/ch07
	$(MAKE) clean --directory=src/ch10
