FC = caf
FCFLAGS = -O3

OBJS = mod_diff.o mod_initial.o mod_parallel.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: tsunami

tsunami: tsunami.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

mod_diff.o: mod_diff.f90
mod_initial.o: mod_initial.f90
mod_parallel.o: mod_parallel.f90

clean:
	$(RM) tsunami *.o *.mod
