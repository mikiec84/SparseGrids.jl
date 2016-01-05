void interp_cc_l(int D, int Q, int nG, int nx,
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
                temp2 *= bf_cc_l( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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

void interp_m_l(int D, int Q, int nG, int nx,
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
                temp2 *= bf_m_l( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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

void interp_m_q(int D, int Q, int nG, int nx,
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
                temp2 *= bf_m_q( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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


void interp_cc_l_big(long D, long Q, long nG, long nx, long nc,
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
                B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1));
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
                // B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1));
                B[q-1+d*(Q+1)] = 1.0-pow((x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1),2);
                // dx = 1.0-pow((mi-1.0)*(x-xij),2);
                // dx = 1.0-((mi-1.0)*fabs(x-xij));
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
