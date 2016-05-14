void w_cc_l(int nG, int D, int Q,
            double *grid, double *A, short *lvl_s, int *lvl_l,
            double *Aold, double *dA, double *w)
{
    double temp, temp2;
    int i, ii, d, q;

    for (q=0;q<Q+1;q++)
    {
        for (i=lvl_l[q]-1;i<=lvl_l[q+1]-2;i++)
            w[i] = A[i] - Aold[i];

        #pragma omp parallel for private(temp,ii,temp2,d)
        for (i=lvl_l[q+1]-1;i<nG;i++)
        {
            temp = 0;
            dA[i] = 0;

            for (ii=lvl_l[q];ii<lvl_l[q+1];ii++)
            {
                temp2 =1;
                for (d=0;d<D;d++)
                {
                    temp2 *= bf_cc_l( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    if (temp2==0.0)
                        break;
                }
                temp += temp2*w[ii-1];
            }
            dA[i] = temp;
        }

        for (i=0;i<nG;i++)
            Aold[i] += dA[i];
    }
}

void w_cc_l_arr(double *grid, int nG, int D, int Q, int nA,
            short *lvl_s, int *lvl_l, double *A,
            double *Aold, double *w)
{
    double temp2;
    int i, ii, iii, d, q;

    for (q=0;q<Q+1;q++)
    {
        for (i=lvl_l[q]-1;i<=lvl_l[q+1]-2;i++)
        {
            for (ii=0;ii<nA;ii++)
                w[i+ii*nG] = A[i+ii*nG] - Aold[i+ii*nG];
        }
        #pragma omp parallel for private(ii,iii,temp2,d)
        for (i=lvl_l[q+1]-1;i<nG;i++)
        {
            for (ii=lvl_l[q];ii<lvl_l[q+1];ii++)
            {
                temp2 =1;
                for (d=0;d<D;d++)
                {
                    temp2 *= bf_cc_l( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    if (temp2==0.0)
                        break;
                }
                for (iii=0;iii<nA;iii++)
                    Aold[i+iii*nG] += temp2*w[ii-1+iii*nG];
            }
        }
    }
}


void interp_l(int D, int Q, int nG, int nx,
                     double * grid, short * lvl_s,
                     double * w,
                     double * xi, double * y,
                     int * nextid)
{
    double temp2;
    int i, ii, d;

    #pragma omp parallel for private(ii,temp2,d)
    for (i=0;i<nx;i++)
    {
        ii=1;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= bf_cc_l( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                if (temp2==0.0)
                    break;
            }
            if (temp2==0.0)
                ii+=1;
            else
            {
                y[i] += temp2*w[ii-1];
                ii=nextid[ii-1];
            }
        }
    }
}


void interp_cc_l_arr(int D, int Q, int nG, int nx, int nA,
                     double * grid, short * lvl_s,
                     double * w,
                     double * xi, double * y,
                     int * nextid)
{
    double temp2;
    int i, ii, iii, d;

    #pragma omp parallel for private(ii,iii,temp2,d)
    for (i=0;i<nx;i++)
    {
        ii=1;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= bf_cc_l( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                if (temp2==0.0)
                    break;
            }
            if (temp2==0.0)
                ii+=1;
            else
            {
                for (iii=0;iii<nA;iii++)
                    y[i+iii*nx] += temp2*w[ii-1+iii*nG];
                ii=nextid[ii-1];
            }
        }
    }
}


void w_cc_q_arr(double *grid, int nG, int D, int Q, int nA,
            short *lvl_s, int *lvl_l, double *A,
            double *Aold, double *w)
{
    double temp2;
    int i, ii, iii, d, q;

    for (q=0;q<Q+1;q++)
    {
        for (i=lvl_l[q]-1;i<=lvl_l[q+1]-2;i++)
        {
            for (ii=0;ii<nA;ii++)
                w[i+ii*nG] = A[i+ii*nG] - Aold[i+ii*nG];
        }
        #pragma omp parallel for private(ii,iii,temp2,d)
        for (i=lvl_l[q+1]-1;i<nG;i++)
        {
            for (ii=lvl_l[q];ii<lvl_l[q+1];ii++)
            {
                temp2 =1;
                for (d=0;d<D;d++)
                {
                    temp2 *= bf_cc_q( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    if (temp2==0.0)
                        break;
                }
                for (iii=0;iii<nA;iii++)
                    Aold[i+iii*nG] += temp2*w[ii-1+iii*nG];
            }
        }
    }
}

void w_cc_q(double *grid, int nG, int D, int Q, short *lvl_s, int *lvl_l, double *A, double *Aold, double *dA, double *w)
{
    double temp, temp2;
    int i, ii, d, q;

    for (q=0;q<Q+1;q++)
    {
        for (i=lvl_l[q]-1;i<=lvl_l[q+1]-2;i++)
            w[i] = A[i] - Aold[i];

        #pragma omp parallel for private(temp,ii,temp2,d)
        for (i=lvl_l[q+1]-1;i<nG;i++)
        {
            temp = 0;
            dA[i] = 0;

            for (ii=lvl_l[q];ii<lvl_l[q+1];ii++)
            {
                temp2 =1;
                for (d=0;d<D;d++)
                {
                    temp2 *= bf_cc_q( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    if (temp2==0.0)
                        break;
                }
                temp += temp2*w[ii-1];
            }
            dA[i] = temp;
        }
        for (i=0;i<nG;i++)
            Aold[i] += dA[i];
    }
}

void interp_cc_q(int D, int Q, int nG, int nx,
                     double * grid, short * lvl_s,
                     double * w,
                     double * xi, double * y,
                     int * nextid)
{
    double temp2;
    int i, ii, d;

    #pragma omp parallel for private(ii,temp2,d)
    for (i=0;i<nx;i++)
    {
        ii=1;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= bf_cc_q( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                if (temp2==0.0)
                    break;
            }
            if (temp2==0.0)
                ii+=1;
            else
            {
                y[i] += temp2*w[ii-1];
                ii=nextid[ii-1];
            }
        }
    }
}

void interp_cc_q_arr(int D, int Q, int nG, int nx, int nA,
                     double * grid, short * lvl_s,
                     double * w,
                     double * xi, double * y,
                     int * nextid)
{
    double temp2;
    int i, ii, iii, d;

    #pragma omp parallel for private(ii,iii,temp2,d)
    for (i=0;i<nx;i++)
    {
        ii=1;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= bf_cc_q( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                if (temp2==0.0)
                    break;
            }

            if (temp2==0.0)
                ii+=1;
            else
            {
                for (iii=0;iii<nA;iii++)
                    y[i+iii*nx] += temp2*w[ii-1+iii*nG];
                ii=nextid[ii-1];
            }
        }
    }
}
