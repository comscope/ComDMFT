include ../../arch.mk

# Implict rules
.SUFFIXES:.f90 .f .o

o_files = \
  flapwmbpt_input_sub.o \

all: 
	$(fortran2python) -m flapwmbpt_input_sub flapwmbpt_input_sub.F

.f90.o:
	$(comp) $(INCLUDE) $<
.f.o:
	$(comp) $(INCLUDE) $<

%.mod: %.f
	$(comp) $<

%.mod: %.f90
	$(comp) $<
clean: 
	rm -rf *.so
	rm -rf *.mod
	rm -rf *.exe
