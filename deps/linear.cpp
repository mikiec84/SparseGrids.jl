inline double linear_bf (double x, double xij, int mi)
{
    double dx;
    if (mi==1)
        return 1.0;
    dx = 1.0-((mi-1.0)*fabs(x-xij));
    if (dx>0)
        return dx;
    else return 0.0;
}


void w_get_l(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
{
    double temp, temp2;
    int i, ii, d, q;

    w[0] = A[0];
    for (i=0;i<nG;i++)
        Aold[i] = A[0];

    for (q=1;q<Q+1;q++)
    {
        for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++)
            w[i] = A[i] - Aold[i];

        #pragma omp parallel for private(temp,ii,temp2,d)
        for (i=lvl_l[q]-1;i<nG;i++)
        {
            temp = 0;
            dA[i] = 0;

            for (ii=lvl_l[q-1];ii<lvl_l[q];ii++)
            {
                temp2 =1;
                for (d=0;d<D;d++)
                {
                    temp2 *= linear_bf( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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


void sparse_interp_l(int D, int Q, int nG, int nx,
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
        ii=2;
        while (ii<=nG)
        {
            temp2 =1.0;
            for (d=0;d<D;d++)
            {
                temp2 *= linear_bf( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
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


void w_get_inv_l(double *grid, int nG, int D, int *lvl_s, int *lvl_l, double *A, int Q, double *Aold, double *dA, double *w)
{
    double temp[nG], temp2;
    int ii, d;
    int q;
    int i,j;

    for (q=0;q<Q+1;q++){
        if (q==0){
            for (i=0;i<=lvl_l[0]-2;i++){
                for (j=0;j<nG;j++){
                    w[i+nG*j] = A[i+nG*j] - Aold[i+nG*j];}}}
        else{
            for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++){
                for (j=0;j<nG;j++){
                    w[i+nG*j] = A[i+nG*j] - Aold[i+nG*j];}}}

        #pragma omp parallel for private(temp,ii,temp2,d,j)
        for (i=0;i<nG;i++)
        {
            for (j=0;j<nG;j++)
                temp[j] = 0;
            if (q == 0)
            {
                for (ii=0;ii<=lvl_l[1]-2;ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= linear_bf( grid[i+d*nG] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
                    for (j=0;j<nG;j++)
                        temp[j] += temp2*w[ii+nG*j];
                }
            }
            else
            {
                for (ii=lvl_l[q-1];ii<lvl_l[q];ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= linear_bf( grid[i+d*nG] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
                    for (j=0;j<nG;j++)
                        temp[j] += temp2*w[ii-1+nG*j];
                }
            }
            for (j=0;j<nG;j++)
                Aold[i+nG*j] += temp[j];
        }
    }
}


void q_get_l(double * xi, int nx, double * grid, int nG, int D, int * lvl_s, int * lvl_l, int Q,double * X)
{
    double temp2;
    int q, i, ii, d;

    for (q=0;q<Q+1;q++)
    {
        #pragma omp parallel for private(ii,temp2,d)
        for (i=0;i<nx;i++)
        {
            if (q == 0)
            {
                for (ii=0;ii<(lvl_l[1]-1);ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= linear_bf( xi[i+d*nx] , grid[ii+d*nG] , lvl_s[ii+d*nG] );

                    X[i+ii*nx] += temp2;
                }
            }
            else
            {
                for (ii=lvl_l[q-1]-1;ii<(lvl_l[q]-1);ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= linear_bf( xi[i+d*nx] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
                    X[i+ii*nx] += temp2;
                }
            }
        }
    }
}

void q_get_l2(double * xi, int nx, double * grid, int nG, int D, int * lvl_s, int * lvl_l, int Q,double * X)
{
    double temp2;
    int q, i, ii, d;

    for (q=0;q<Q+1;q++)
    {
        #pragma omp parallel for private(ii,temp2,d)
        for (i=0;i<nx;i++)
        {
            if (q == 0)
            {
                X[i] += 1.0;
            }
            else
            {
                for (ii=lvl_l[q-1]-1;ii<(lvl_l[q]-1);ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                    {
                        temp2 *= linear_bf( xi[i+d*nx] , grid[ii+d*nG] , lvl_s[ii+d*nG] );
                        if (temp2==0.0)
                            break;
                    }
                    X[i+ii*nx] += temp2;
                }
            }
        }
    }
}



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
