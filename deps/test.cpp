void interp_cc_l_big1(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, int* combl, int* combu)
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
                B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i +d*nx]-dxi_cc(q,j))*(mi-1));
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)];
            for (d=D-2;d>0;d--)
            {
                b*=B[combs[ii+d*nc]-1+d*(Q+1)];
                id1 = id1*dM_cc(combs[ii+d*nc]) + J[combs[ii+(d)*nc]-1+(d)*(Q+1)];
            }
            id1 = id1*dM_cc(combs[ii])+(J[combs[ii]-1])+combl[ii]-1;
            y[i]+=b*w[id1];
        }
    }
}
