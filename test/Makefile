sparsemake: funcs.cpp
	gcc -Wall -O3 -fopenmp -lgomp -fpic -c funcs.cpp
	gcc -shared -O3 -fopenmp -lgomp -Wl,-soname,libsparse1.so -o libsparse1.so   funcs.o
	sudo cp libsparse1.so /usr/lib/libsparse1.so
	sudo rm funcs.o
	nm /usr/lib/libsparse1.so |grep interp
