GCC=gcc-14

export GCC

all:
	cd ./SEmain/ && make
	cd ./SEmain/ && make install
	cd ./src && make

install: 
	cd ./src && make install
	
clean:
	cd ./SEmain/SEbasic/src && make clean
	cd ./SEmain/SEmath/src && make clean
	cd ./SEmain/SEutils/src && make clean
	rm -f libSE.a test
	cd ./src && make clean
	rm -f ../lib/libSE.a
	
