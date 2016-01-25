
void interp_cc_l_big_s(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, int* combl, int* combu, long* hGji)
{
    long i,ii,iii,j,d,q,id1;
    long J[D*(Q+1)];
    long long hid;
    double B[D*(Q+1)],b;

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


void interp_cc_l_big1_s(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, short* combs_dM, int* combl, int* combu, long* hGji)
{
    long i,ii,j,d,q,id1;
    long J[D*(Q+1)];
    double B[D*(Q+1)],b;

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
                J[q-1+d*(Q+1)] = j ;
                B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1));
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]*B[combs[ii]-1];
            id1 = (J[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]-1);
            for (d=D-2;d>0;d--)
            {
                b*=B[combs[ii+d*nc]-1+d*(Q+1)];
                id1 = id1*combs_dM[ii+d*nc] + (J[combs[ii+(d)*nc]-1+(d)*(Q+1)]-1);
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1]-1)+combl[ii]-1;
            y[i]+=b*w[id1];
        }
    }
}


void interp_cc_l_big1_p(long D, long Q, long nG, long nx, long nc,
    double* x, double* w, double* y,
    short* combs, short* combs_dM, int* combl, int* combu, long* hGji)
{
    long i,ii,j,d,q,id1;
    long J[D*(Q+1)];
    double B[D*(Q+1)],b;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1)
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
                J[q-1+d*(Q+1)] = j ;
                B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i +d*nx]-dxi_cc(q,j))*(M_cc(q)-1));
            }
        }

        for (ii=0;ii<nc;ii++)
        {
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]*B[combs[ii]-1];
            id1 = (J[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)]-1);
            for (d=D-2;d>0;d--)
            {
                b*=B[combs[ii+d*nc]-1+d*(Q+1)];
                id1 = id1*combs_dM[ii+d*nc] + (J[combs[ii+(d)*nc]-1+(d)*(Q+1)]-1);
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1]-1)+combl[ii]-1;
            y[i]+=b*w[id1];
        }
    }
}
