FC = gfortran
FCFLAGS = -O3

OBJS = mod_diff.o mod_initial.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: tsunami

tsunami: tsunami.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

clean:
	$(RM) tsunami *.o *.mod
