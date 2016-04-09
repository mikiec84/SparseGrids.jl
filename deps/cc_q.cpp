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


void interp_cc_q_big(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, short* combs_dM, int* combl, int* combu)
{
    long i,ii,j,d,q,id1,dmi;
    long J[D][Q+1];
    double B[D][Q+1],b;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi)
    for (i=0;i<nx;i++)
    {
        for (d=0;d<D;d++)
        {
            J[d][0] = 0;
            B[d][0] = 1.0;
            J[d][1] = (x[i+d*nx]>0.5);
            B[d][1] = 1.0-pow((x[i+d*nx]-J[d][1])*2,2);
        }
        for (q=2;q<Q+1;q++)
        {
            dmi = pow(2,q-1);
            for (d=0;d<D;d++)
            {
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[d][q] = j-1;
                B[d][q] = 1.0-pow(2*(x[i+d*nx]*dmi-(j-0.5)),2);
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[D-1][combs[ii+(D-1)*nc]-1] * B[0][combs[ii]-1];
            id1 = J[D-1][combs[ii+(D-1)*nc]-1];
            for (d=D-2;d>0;d--)
            {
                b*=B[d][combs[ii+d*nc]-1];
                id1 = id1*combs_dM[ii+d*nc] + J[d][combs[ii+(d)*nc]-1];
            }
            id1 = id1*combs_dM[ii]+(J[0][combs[ii]-1])+combl[ii]-1;
            y[i]+=b*w[id1];
        }
    }
}


void interp_cc_q_big_arr(long D, long Q, long nG, long nx, long nc, long nA,
    double* x, double* w, double* y,
    short* combs, short* combs_dM, int* combl, int* combu)
{
    long i,ii,j,d,q,id1,dmi;
    long J[D][Q+1];
    double B[D][Q+1],b;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi)
    for (i=0;i<nx;i++)
    {
        for (d=0;d<D;d++)
        {
            J[d][0] = 0;
            B[d][0] = 1.0;
            J[d][1] = (x[i+d*nx]>0.5);
            B[d][1] = 1.0-pow((x[i+d*nx]-J[d][1])*2,2);
        }
        for (q=2;q<Q+1;q++)
        {
            dmi = pow(2,q-1);
            for (d=0;d<D;d++)
            {
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[d][q] = j-1;
                B[d][q] = 1.0-pow(2*(x[i+d*nx]*dmi-(j-0.5)),2);
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[D-1][combs[ii+(D-1)*nc]-1] * B[0][combs[ii]-1];
            id1 = J[D-1][combs[ii+(D-1)*nc]-1];
            for (d=D-2;d>0;d--)
            {
                b*=B[d][combs[ii+d*nc]-1];
                id1 = id1*combs_dM[ii+d*nc] + J[d][combs[ii+(d)*nc]-1];
            }
            id1 = id1*combs_dM[ii]+(J[0][combs[ii]-1])+combl[ii]-1;
            for (d=0;d<nA;d++)
                y[i+d*nx] += b*w[id1+d*nG];
        }
    }
}
