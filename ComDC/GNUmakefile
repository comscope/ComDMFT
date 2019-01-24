include ../arch.mk
comp = $(PF90) $(FPPFLAGS) -c $(compfl)
link = $(PF90)  $(compfl)
LIB = $(LAPACK_LIB)

# Implict rules
.SUFFIXES:.mod .F .o

o_files = \
  comdc.o \
  comdc_mod.o \

mod_files = \
  comdc_mod.mod

TARGET = ComDC

$(TARGET): $(mod_files) $(o_files) 
	$(link) $(o_files) -o $(TARGET) $(LIB)
	cp $(TARGET) ../bin
.f90.o:
	$(comp) $(INCLUDE) $<
.F.o:
	$(comp) $(INCLUDE) $<
%.mod: %.F
	$(comp) $<
%.mod: %.f90
	$(comp) $<
clean: 
	rm -rf *.o 
	rm -rf *.mod
	rm -rf $(TARGET)

