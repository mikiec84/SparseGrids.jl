sparsemake: funcs.cpp
	gcc -Wall -O3 -fopenmp -lgomp -fpic -c funcs.cpp
	gcc -shared -O3 -fopenmp -lgomp -Wl,-soname,libsparse.so -o libsparse.so   funcs.o
	sudo cp libsparse.so /usr/lib/libsparse.so
	sudo rm funcs.o
