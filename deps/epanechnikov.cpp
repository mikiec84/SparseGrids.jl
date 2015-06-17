

inline double epanechnikov (double x, double xij, int mi)
{
    if (mi==1)
        return 1;
    else if (abs(x-xij)<(1./(mi-1)))
        return (1-pow(((mi)-1)*(x-xij),2));
    else return 0;
}




void w_get_epa(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
{
    double temp, temp2;
    int ii, d;
    int q;
    int i;

    for (q=0;q<Q+1;q++){
        if (q==0)
            w[0]=A[0];
        else{
            for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++){
                w[i] = A[i] - Aold[i];
            }
        }

        #pragma omp parallel for private(temp,ii,temp2,d)
        for (i=0;i<nG;i++){
            temp = 0;
            dA[i] = 0;
            if (q == 0){
                temp = w[0];
            }
            else{
                for (ii=lvl_l[q-1];ii<lvl_l[q];ii++){
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= epanechnikov( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    temp += temp2*w[ii-1];
                }
            }
            dA[i] = temp;
        }

        for (i=0;i<nG;i++)
            Aold[i] += dA[i];
    }

}


void sparse_interp_epa(double * xi, int nx, double * grid, int nG, int D,
    int * lvl_s, int * lvl_l, double * A, int Q, double * w,
    double * xold, double * dx)
{
    double temp, temp2;
    int q, i, ii, d;



    for (q=0;q<Q+1;q++)
    {
        #pragma omp parallel for private(temp,ii,temp2,d)
        for (i=0;i<nx;i++)
        {
            temp = 0;
            dx[i] = 0;
            if (q == 0)
            {
                temp = w[0];
            }
            else
            {
                for (ii=lvl_l[q-1];ii<(lvl_l[q]);ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                    {
                        temp2 *= epanechnikov( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    }
                    temp += temp2*w[ii-1];
                }
            }
            dx[i] = temp;
        }

        #pragma omp parallel for
        for (i=0;i<nx;i++)
            xold[i] += dx[i];

    }
}
