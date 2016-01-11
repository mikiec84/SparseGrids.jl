sparsemake: libsparse.cpp
	gcc -Wall -O3  -fopenmp -lgomp -fpic -c libsparse.cpp
	gcc -shared -O3 -fopenmp -lgomp -Wl,-soname,libsparse.so -o libsparse.so   libsparse.o
	sudo cp libsparse.so /usr/lib/libsparse.so
	sudo rm libsparse.o
