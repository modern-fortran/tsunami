# tsunami Makefile

.PHONY: all clean

all: tsunami

tsunami: tsunami.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	$(RM) tsunami
