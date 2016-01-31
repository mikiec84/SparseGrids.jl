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
        // y[i] += w[0];
        // ii=2;
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


void interp_cc_q_big(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, short* combs_dM, int* combl, int* combu)
{
    long i,ii,j,d,q,id1,dmi,mi;
    long J[D*(Q+1)];
    double B[D*(Q+1)],b;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi)
    for (i=0;i<nx;i++)
    {
        for (q=1;q<=Q+1;q++)
        {
            dmi = dM_cc(q);
            mi = M_cc(q);
            for (d=0;d<D;d++)
            {
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[q-1+d*(Q+1)] = j-1;
                B[q-1+d*(Q+1)] = 1.0-pow((x[i +d*nx]-dxi_cc(q,j))*(mi-1),2);
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)];
            for (d=D-2;d>0;d--)
            {
                b*=B[combs[ii+d*nc]-1+d*(Q+1)];
                id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(Q+1)];
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1])+combl[ii]-1;
            y[i]+=b*w[id1];
        }
    }
}


void interp_cc_q_big_arr(long D, long Q, long nG, long nx, long nc, long nA,
    double* x, double* w, double* y,
    short* combs, short* combs_dM, int* combl, int* combu)
{
    long i,ii,j,d,q,id1,dmi,mi;
    long J[D*(Q+1)];
    double B[D*(Q+1)],b;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi)
    for (i=0;i<nx;i++)
    {
        for (q=1;q<=Q+1;q++)
        {
            dmi = dM_cc(q);
            mi = M_cc(q);
            for (d=0;d<D;d++)
            {
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[q-1+d*(Q+1)] = j-1;
                B[q-1+d*(Q+1)] = 1.0-pow((x[i +d*nx]-dxi_cc(q,j))*(mi-1),2);
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)];
            for (d=D-2;d>0;d--)
            {
                b*=B[combs[ii+d*nc]-1+d*(Q+1)];
                id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(Q+1)];
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1])+combl[ii]-1;
            for (d=0;d<nA;d++)
                y[i+d*nx] += b*w[id1+d*nG];
        }
    }
}
