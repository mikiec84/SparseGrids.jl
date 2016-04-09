sparsemake: libsparse.cpp
	gcc -Wall -O3  -fopenmp -lgomp -fpic -g -c libsparse.cpp
	gcc -shared -O3 -fopenmp -g -lgomp -Wl,-soname,libsparse.so -o libsparse.so   libsparse.o
