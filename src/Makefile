IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -O3 -Wunused-variable

ODIR=obj

LIBS=-lm

_DEPS = functions.h strucs.h vec.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

SOURCES = *.c

_OBJ = fittness.o initiate_genotypes.o main.o micutc.o mutate.o various.o \
       copy_gene.o mate.o setup_mm.o allocate_memory.o read_topology.o write_mme.o \
       write_top.o get_parameters.o read_qme.o read_coords.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fogal:	$(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY:	clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ fogal

install: fogal
	cp fogal ~/bin/
