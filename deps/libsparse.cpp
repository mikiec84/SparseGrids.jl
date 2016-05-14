#include <cmath>
#include <omp.h>
#include <thread>
#include <iostream>
#include <emmintrin.h>
#include <stdint.h>
#include <cassert>
#include <nmmintrin.h>
#include <tmmintrin.h>
using namespace std;


#include "griddefs.cpp"
#include "intersect.cpp"
#include "cc_l.cpp"
#include "cc_q.cpp"
#include "test.cpp"


void c_updateU(void *Fval, void *Jval, void *X, void *XP, void *SP, void *ProbWeights, long N,void(*f)(void *z, void *a, void *b, void *c, void *d, void *e, const long i))
{
    int j;
    #pragma omp parallel for
    for (j=1;j<=N;j++)
        f(Fval,Jval,X,XP,SP,ProbWeights,j);
}
