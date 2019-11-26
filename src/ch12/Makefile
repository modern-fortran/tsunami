# tsunami Makefile

FC = caf
FCFLAGS = -O3 -cpp

.PHONY: all clean
.SUFFIXES: .f90 .o

OBJS = mod_diff.o mod_field.o mod_io.o mod_kinds.o mod_parallel.o

all: tsunami

tsunami: tsunami.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

mod_diff.o: mod_diff.f90 mod_kinds.o 
mod_field.o: mod_field.f90 mod_io.o mod_kinds.o mod_parallel.o
mod_kinds.o: mod_kinds.f90
mod_parallel.o: mod_parallel.f90 mod_kinds.o

clean:
	$(RM) tsunami *.o *.mod
