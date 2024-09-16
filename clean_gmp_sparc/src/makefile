CC=gcc
CFLAGS= -g -Iinclude

OBJSC = gmp_deriv_descriptor.o helper.o local_type_tools.o solid_harmonics.o surface_harmonics.o TEST_GMP_SPARC.o tools_mlff.o
LDLIBS = -lm

LIBBASE = ../lib/gmp

all: GMP

GMP: $(OBJSC)
	$(CC) $(CFLAGS) -o $(LIBBASE) $^ $(LDLIBS)

.PHONY: clean

clean:
	rm -f $(OBJSC) $(LIBBASE)
