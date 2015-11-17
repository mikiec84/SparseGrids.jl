inline double nbfirstXi(int mi)
{
    return double(1)/double(mi+1);
}

inline double nblastXi(int mi)
{
    return double(mi)/double(mi+1);
}

double basis_func_nb (double x, double xij, int mi)
{
    if (mi==1)
        return 1.0;
    else
    {
        if (xij==nbfirstXi(mi))
        {
            if (x<2./(mi+1.))
                return 2.-(mi+1.)*x;
            else
                return 0.0;
        }
        if (xij==nblastXi(mi))
        {
            if (x>(mi-1.)/(mi+1.))
                return (mi+1.)*x -mi + 1.;
            else
                return 0.0;
        }
        else
        {
            if (fabs(x-xij)<1/(mi+1.))
                return 1.-(mi+1.)*fabs(x-xij);
            else
                return 0.0;

        }

    }
}

void w_get_nb(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
{
    double temp, temp2;
    int ii, d;
    int q;
    int i;

    for (q=0;q<Q+1;q++){
        if (q==0)
        {
                for (i=0;i<=lvl_l[0]-2;i++)
                {
                    w[i] = A[i] - Aold[i];
                }
        }
        else
        {
            for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++)
            {
                w[i] = A[i] - Aold[i];
            }
        }

        #pragma omp parallel for private(temp,ii,temp2,d)
        for (i=0;i<nG;i++)
        {
            temp = 0;
            dA[i] = 0;
            if (q == 0)
            {
                for (ii=0;ii<=lvl_l[1]-2;ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= basis_func_nb( grid[i+d*nG] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
                    temp += temp2*w[ii];
                }
            }
            else{
                for (ii=lvl_l[q-1];ii<lvl_l[q];ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= basis_func_nb( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    temp += temp2*w[ii-1];
                }
            }
            dA[i] = temp;
        }
        for (i=0;i<nG;i++)
            Aold[i] += dA[i];
    }
}

void w_get_inv_nb(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
{
    double temp[nG], temp2;
    int ii, d;
    int q;
    int i,j;

    for (q=0;q<Q+1;q++){
        if (q==0)
        {
            for (i=0;i<=lvl_l[0]-2;i++)
            {
                for (j=0;j<nG;j++)
                {
                    w[i+nG*j] = A[i+nG*j] - Aold[i+nG*j];
                }

            }
        }
        else
        {
            for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++)
            {
                for (j=0;j<nG;j++)
                {
                    w[i+nG*j] = A[i+nG*j] - Aold[i+nG*j];
                }
            }
        }

        #pragma omp parallel for private(temp,ii,temp2,d,j)
        for (i=0;i<nG;i++)
        {
            for (j=0;j<nG;j++)
            {
                temp[j] = 0;
            }


            if (q == 0)
            {
                for (ii=0;ii<=lvl_l[1]-2;ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= basis_func_nb( grid[i+d*nG] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
                    for (j=0;j<nG;j++)
                    {
                        temp[j] += temp2*w[ii+nG*j];
                    }
                }
            }
            else
            {
                for (ii=lvl_l[q-1];ii<lvl_l[q];ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= basis_func_nb( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    for (j=0;j<nG;j++)
                    {
                        temp[j] += temp2*w[ii-1+nG*j];
                    }
                }
            }
            for (j=0;j<nG;j++)
                Aold[i+nG*j] += temp[j];
        }
    }
}


// void w_get_nb(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
// {
//     double temp, temp2;
//     int ii, d;
//     int q;
//     int i;

//     for (q=0;q<Q+1;q++){
//         if (q==0)
//             w[0]=A[0];
//         else{
//             for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++){
//                 w[i] = A[i] - Aold[i];
//             }
//         }

//         #pragma omp parallel for private(temp,ii,temp2,d)
//         for (i=0;i<nG;i++){
//             temp = 0;
//             dA[i] = 0;
//             if (q == 0){
//                 temp = w[0];
//             }
//             else{
//                 for (ii=lvl_l[q-1];ii<lvl_l[q];ii++){
//                     temp2 =1;
//                     for (d=0;d<D;d++)
//                         temp2 *= basis_func_nb( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
//                     temp += temp2*w[ii-1];
//                 }
//             }
//             dA[i] = temp;
//         }

//         for (i=0;i<nG;i++)
//             Aold[i] += dA[i];
//     }

// }

void sparse_interp_nb(int D, int Q, int nG, int nx,
    double * grid, int * lvl_s, int * lvl_l,
    double * A,  double * w,
    double * xi, double * xold, double * dx)
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
                        temp2 *= basis_func_nb( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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
