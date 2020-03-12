FCOMPILER=gfortran -ff2c -I. -ffixed-line-length-none -fno-automatic -fno-range-check -Wall -g
FLINK= gfortran -ff2c -mcmodel=large

BINDIR=../bin
OBJDIR=../obj
INCDIR=../include
LIBDIR=/usr/local/lib

LIBS=$(LIBDIR)/libumlm.a $(LIBDIR)/libgeo.a $(LIBDIR)/libio.a

OBJ=$(OBJDIR)/main.o $(OBJDIR)/ylmpav.o $(OBJDIR)/forward_modelling.o

$(OBJDIR)/%.o: %.f
	$(FCOMPILER) -c -I$(INCDIR) $^; mv *.o $(OBJDIR)

main: $(OBJ) 
	$(FLINK) $^ $(LIBS) -o $(BINDIR)/$@