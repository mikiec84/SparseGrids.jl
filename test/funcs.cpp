#include <cmath>
#include <omp.h>
#include <iostream>
#include <map>
using namespace std;



int M(long i)
{
    if (i==1)
        return 1;
    else
        return pow(2,i-1)+1;
}

int dM(long i)
{
    if (i==1)
        return 1;
    else if (i==2)
        return 2;
    else
        return pow(2,i-2);
}

double dxi(long i, long j)
{
    if (i==1)
        return 0.5;
    else if (i==2)
        if (j==1)
            return 0.0;
        else
            return 1.0;
    else
        return (j*2.0-1.0)/(M(i)-1);
}

int indexshift(long i, long j, long q)
{
    if (i==1)
        return (M(q)-1)/2+1;
    else if (i==2)
        if (j==1)
            return 1;
        else
            return M(q);
    else
        return ((M(q)-1)*(j*2-1))/(M(i)-1)+1;

}



void interp_l_big(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y, long* combs, long* combl, long* combu, long* hGji)
{
    long i,ii,iii,j,d,q,id1;
    long J[D*(Q+1)];
    long long hid;
    double B[D*(Q+1)],b;


    #pragma omp parallel for private(d,q,j,J,B,ii,b,hid,id1,iii)
    for (i=0;i<nx;i++)
    {
        for (d=0;d<D;d++)
        {
            for (q=1;q<=Q+1;q++)
            {
                j = round(x[i+d*nx]*dM(q)+0.5);
                if (j<1)
                    j=1;
                else if (j>dM(q))
                    j = dM(q);
                J[q-1+d*(Q+1)] = indexshift(q,j,Q+1);
                B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i +d*nx]-dxi(q,j))*(M(q)-1));

            }
        }


        for (ii=0;ii<nc;ii++)
        {
            b=1.0;
            hid = 0;
            for (d=0;d<D;d++)
            {
                b*=B[combs[ii+d*nc]-1+d*(Q+1)];
                hid += J[combs[ii+d*nc]-1+d*(Q+1)];
                hid *=17;
            }
            id1 = 0;
            for (iii=combl[ii]-1;iii<combu[ii];iii++)
            {
                if (hGji[iii]==hid)
                {
                    id1 = iii;
                    break;
                }
            }
            y[i]+=b*w[id1];
        }
    }
}
