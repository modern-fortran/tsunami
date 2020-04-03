FC = caf
FCFLAGS = -O3

OBJS = mod_diff.o mod_field.o mod_io.o mod_parallel.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: tsunami

tsunami: tsunami.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

mod_field.o: mod_field.f90 mod_io.o mod_parallel.o

clean:
	$(RM) tsunami *.o *.mod
