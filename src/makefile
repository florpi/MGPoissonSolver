IDIR =../include
CC=g++
CFLAGS=-I$(IDIR)

ODIR=obj

_DEPS=main.h grid.h multigrid.h particle.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =main.o grid.o multigrid.o particle.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -O2 -o $@ $< $(CFLAGS)
acceleration: $(OBJ)
	g++  -o $@ $^ $(CFLAGS) 

.PHONY:	clean
clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
