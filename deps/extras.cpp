
int ipow(int x, int n)
{
    int out;
    out = x;
    for (int i = 1;i<n;i++)
        out*=x;
    return out;
}

int ij_Mi(int i)
{
    if (i==1)
        return 1;
    else
        return ipow(2,i-1)+1;
}


double ij_xi(int i, int j)
{
    if (i==1)
        return 0.5;
    else
        return (double(j-1))/double(ij_Mi(i)-1);
}


double ij_linear_bf (double x, int i, int j)
{
    if (i==1)
        return 1.0;
    int mi;
    mi = ij_Mi(i);
    double xij=ij_xi(i,j);
    if (abs(x-xij)<(1./(mi-1)))
        return (1-((mi)-1)*abs(x-xij));
    else return 0.0;
}

void w_get_l_ij(int Q, int nG, int D,double *grid,  int *lvl_i, int *lvl_j, int *lvl_l, double *A  , double *w)
{
    double temp, temp2;
    int ii, d;
    int q;
    int i;
    double *dA,*Aold;
    dA = new double [nG];
    Aold = new double [nG];

    for (q=0;q<Q+1;q++){
        if (q==0)
                for (i=0;i<=lvl_l[0]-2;i++)
                    w[i] = A[i] - Aold[i];
        else
            for (i=lvl_l[q-1]-1;i<=lvl_l[q]-2;i++)
                w[i] = A[i] - Aold[i];

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
                        temp2 *= linear_bf( grid[i+d*nG] , lvl_i[ii+d*nG] , lvl_j[ii+d*nG] );
                    temp += temp2*w[ii];
                }
            }
            else{
                for (ii=lvl_l[q-1];ii<lvl_l[q];ii++)
                {
                    temp2 =1;
                    for (d=0;d<D;d++)
                        temp2 *= linear_bf( grid[i+d*nG] , lvl_i[ii-1+d*nG] , lvl_j[ii-1+d*nG] );
                    temp += temp2*w[ii-1];
                }
            }
            dA[i] = temp;
        }
        for (i=0;i<nG;i++)
            Aold[i] += dA[i];
    }
    delete[] dA;
    delete[] Aold;
}


