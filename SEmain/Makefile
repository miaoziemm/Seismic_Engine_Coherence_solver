
DIR_BASIC=./SEbasic/src
DIR_MATH=./SEmath/src
DIR_UTILS=./SEutils/src

all:
	cd $(DIR_BASIC) && make
	cd $(DIR_MATH) && make
	cd $(DIR_UTILS) && make
	ar -r libSE.a $(DIR_BASIC)/*.o $(DIR_MATH)/*.o $(DIR_UTILS)/*.o

test: test.c
	gcc -o test test.c -L. -lSE -lm

install: all
	cp libSE.a ../lib

clean:
	cd $(DIR_BASIC) && make clean
	cd $(DIR_MATH) && make clean
	cd $(DIR_UTILS) && make clean
	rm -f libSE.a test



