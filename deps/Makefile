sparsemake: libsparse.cpp
	gcc -Wall -std=c++11 -mavx -O3  -pthread -fopenmp -lgomp -fpic -g -c libsparse.cpp
	gcc -shared -O3 -pthread -fopenmp -g -lgomp -Wl,-soname,libsparse.so -o libsparse.so   libsparse.o
