# GCC=gcc

INTERFLAGS=-lm

INCLUDES=-I../include/ -I../../SEbasic/include/

C_LIB=SEmath.c SEmaths.c SEsufft.c SEsort.c SEv3d_m3d.c SEv4d_m4d.c
OBJ=SEmath.o SEmaths.o SEsufft.o SEsort.o SEv3d_m3d.o SEv4d_m4d.o

$(OBJ): $(C_LIB)
	$(GCC) -c $(C_LIB) $(INCLUDES) $(INTERFLAGS)

clean:
	rm -f *.o










