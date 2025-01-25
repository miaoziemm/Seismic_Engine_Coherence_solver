
FLAGS=-lm -fopenmp -O3

LIB=-L../lib -lSE

INCLUDES=-I../SEmain/include 

CFILE=SE_coherence_solver.c SE_coherence.c SE_coherence_math.c SE_coherence_dip.c SE_coherence_func.c
OBJ=SE_coherence_solver.o SE_coherence.o SE_coherence_math.o SE_coherence_dip.o SE_coherence_func.o

APP=se_coherence_solver se_coherence_cut_cube se_coherence_cut_slice se_layer_generator se_interplatation_2d 
 
all: $(APP)



se_coherence_cut_cube: $(OBJ) SE_coherence_cut_cube.c
	$(GCC) -o se_coherence_cut_cube SE_coherence_cut_cube.c $(FLAGS) $(LIB)

se_coherence_cut_slice: $(OBJ) SE_coherence_cut_slice.c
	$(GCC) -o se_coherence_cut_slice SE_coherence_cut_slice.c $(FLAGS) $(LIB)

se_layer_generator: $(OBJ) SE_layer_generator.c
	$(GCC) -o se_layer_generator SE_layer_generator.c $(FLAGS) $(LIB)

se_interplatation_2d: $(OBJ) SE_interplatation_2d.c
	$(GCC) -o se_interplatation_2d SE_interplatation_2d.c $(FLAGS) $(LIB)

se_coherence_solver: $(OBJ)
	$(GCC) -o se_coherence_solver $(OBJ) $(FLAGS) $(LIB)
	
$(OBJ): $(CFILE)
	$(GCC) -c $(CFILE) $(FLAGS) $(INCLUDES)

help:
	cd help; make

install: all
	cp $(APP) ../bin

clean:
	rm -f *.o $(APP)
