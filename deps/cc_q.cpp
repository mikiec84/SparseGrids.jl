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
            y[i] += temp2*w[ii-1];
            if (temp2==0.0)
                ii+=1;
            else
                ii=nextid[ii-1];
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
            for (iii=0;iii<nA;iii++)
                y[i+iii*nx] += temp2*w[ii-1+iii*nG];
            if (temp2==0.0)
                ii+=1;
            else
                ii=nextid[ii-1];
        }
    }
}


void interp_cc_q_big(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, int* combl, int* combu, long* hGji)
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
                j = round(x[i+d*nx]*dM_cc(q)+0.5);
                if (j<1)
                    j=1;
                else if (j>dM_cc(q))
                    j = dM_cc(q);
                J[q-1+d*(Q+1)] = indexshift_cc(q,j,Q+1);
                B[q-1+d*(Q+1)] = 1.0-pow((x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1),2);
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
                hid *= 17;
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



void interp_cc_q_big_arr(long D, long Q, long nG, long nx, long nc, long nA,
    double* x, double* w, double* y,
    short* combs, int* combl, int* combu, long* hGji)
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
                j = round(x[i+d*nx]*dM_cc(q)+0.5);
                if (j<1)
                    j=1;
                else if (j>dM_cc(q))
                    j = dM_cc(q);
                J[q-1+d*(Q+1)] = indexshift_cc(q,j,Q+1);
                B[q-1+d*(Q+1)] = 1.0-pow((x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1),2);
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

            // y[i]+=b*w[id1];
            for (iii=0;iii<nA;iii++)
                y[i+iii*nx] += b*w[id1+iii*nG];
        }
    }
}
