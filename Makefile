# tsunami Makefile

.PHONY: all clean

all: ch02 ch03 ch04 ch07 ch08 ch10 ch12 final

ch02:
	$(MAKE) --directory=src/$@

ch03:
	$(MAKE) --directory=src/$@

ch04:
	$(MAKE) --directory=src/$@

ch07:
	$(MAKE) --directory=src/$@

ch08:
	$(MAKE) --directory=src/$@

ch10:
	$(MAKE) --directory=src/$@

ch12:
	$(MAKE) --directory=src/$@

final:
	$(MAKE) --directory=src/$@

clean:
	$(MAKE) clean --directory=src/ch02
	$(MAKE) clean --directory=src/ch03
	$(MAKE) clean --directory=src/ch04
	$(MAKE) clean --directory=src/ch07
	$(MAKE) clean --directory=src/ch08
	$(MAKE) clean --directory=src/ch10
	$(MAKE) clean --directory=src/ch12
	$(MAKE) clean --directory=src/final
