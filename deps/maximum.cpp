double basis_func_max (double x, double xij, int mi)
{
    if (fabs(x-xij)<1/(mi-1.0))
        return 1.-(mi-1.)*fabs(x-xij);
    else
        return 0.0;
}

void w_get_max(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
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
                    {
                        temp2 *= basis_func_max( grid[i+d*nG] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
                        if (temp2==0.0)
                            break;
                    }
                    temp += temp2*w[ii];
                }
            }
            else{
                for (ii=lvl_l[q-1];ii<lvl_l[q];ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                    {
                        temp2 *= basis_func_max( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                        if (temp2==0.0)
                            break;
                    }
                    temp += temp2*w[ii-1];
                }
            }
            dA[i] = temp;
        }
        for (i=0;i<nG;i++)
            Aold[i] += dA[i];
    }
}

void w_get_inv_max(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
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
                        temp2 *= basis_func_max( grid[i+d*nG] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
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
                        temp2 *= basis_func_max( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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


void sparse_interp_max(int D, int Q, int nG, int nx,
                      double * grid, int * lvl_s, int * lvl_l,
                      double * A,  double * w,
                      double * xi, double * y,
                      int * nextid)
{
    double temp2;
    int i, ii, d;

    #pragma omp parallel for private(ii,temp2,d)
    for (i=0;i<nx;i++)
    {
        y[i] = 0;
        ii=1;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= basis_func_max( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                if (temp2==0.0)
                    break;
            }
            y[i] += temp2*w[ii-1];
            ii++;
        }
    }
}



void sparse_interp_max2(int D, int Q, int nG, int nx,
                      double * grid, int * lvl_s, int * lvl_l,
                      double * A,  double * w,
                      double * xi, double * y,
                      int * nextid)
{
    double temp2;
    int i, ii, d;

    #pragma omp parallel for private(ii,temp2,d)
    for (i=0;i<nx;i++)
    {
        y[i] += w[0];
        ii=1;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= basis_func_max( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                if (temp2==0.0)
                    break;
            }
            y[i] += temp2*w[ii-1];
            if (temp2==0.0)
                ii+=1;
            else
                ii=nextid[ii-1];
        }
    }
}
