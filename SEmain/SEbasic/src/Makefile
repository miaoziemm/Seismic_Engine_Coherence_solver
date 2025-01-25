# GCC=gcc

INTERFLAGS=-lm

INCLUDES=-I ../include/

C_LIB_basic=SEbasic.c SEprint.c SEcomplex.c SEsys.c
OBJ_basic=SEbasic.o SEprint.o SEcomplex.o SEsys.o

$(OBJ_basic): $(C_LIB_basic)
	$(GCC) -c $(C_LIB_basic) $(INCLUDES) $(INTERFLAGS)

clean:
	rm -f *.o










