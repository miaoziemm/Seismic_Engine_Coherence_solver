# GCC=gcc

INTERFLAGS=-lm

INCLUDES=-I../include/ -I../../SEbasic/include/ -I../../SEmath/include/

C_LIB=SE_interpolation.c SE_agc.c SE_bandpass.c SE_vel_cfg.c
OBJ=SE_interpolation.o SE_agc.o SE_bandpass.o SE_vel_cfg.o

$(OBJ): $(C_LIB)
	$(GCC) -c $(C_LIB) $(INCLUDES) $(INTERFLAGS)

clean:
	rm -f *.o










