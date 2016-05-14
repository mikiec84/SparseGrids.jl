void interp_l(long D, long L, long nx, long nc, double* x, double* w, double* Y, short* combs, short* combs_dM, int* covers_loc)
{
    long i,ii,j,d,q,id1,dmi,mi;
    long J[D*(L+1)];
    double B[D*(L+1)],b,y;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi,y)
    for (i=0;i<nx;i++){
        y = 0.0;
        for (d=0;d<D;d++){
            J[d*(L+1)] = 0;
            J[d*(L+1)+1] = (x[i+d*nx]>0.5);
            B[d*(L+1)] = 1.0;
            B[d*(L+1)+1] = fabs(1.0-fabs(x[i+d*nx]-J[d*(L+1)+1])*2);
        }
        for (q=3;q<=L+1;q++){
            dmi = pow(2,q-2);
            mi = dmi*2+1;
            for (d=0;d<D;d++){
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[q-1+d*(L+1)] = j-1;
                B[q-1+d*(L+1)] = fabs(1.0-fabs(x[i+d*nx]-(j*2.0-1.0)/(mi-1))*(mi-1));
            }
        }
        for (ii=0;ii<nc;ii++){
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)];
            for (d=D-2;d>0;d--){
                b*=B[combs[ii+d*nc]-1+d*(L+1)];
                id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(L+1)];
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1])+covers_loc[ii]-1;
            y+=b*w[id1];
        }
        Y[i] = y;
    }
}

void interp_l(long D, long L, long nG, long nx, long nc, long nA, double* x, double* w, double* y, short* combs, short* combs_dM, int* covers_loc)
{
    long i,ii,j,d,q,id1,dmi,mi;
    long J[D*(L+1)];
    double B[D*(L+1)],b;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi)
    for (i=0;i<nx;i++){
        for (d=0;d<D;d++){
            J[d*(L+1)] = 0;
            J[d*(L+1)+1] = (x[i+d*nx]>0.5);
            B[d*(L+1)] = 1.0;
            B[d*(L+1)+1] = fabs(1.0-fabs(x[i+d*nx]-J[d*(L+1)+1])*2);
        }
        for (q=3;q<=L+1;q++){
            dmi = pow(2,q-2);
            mi = dmi*2+1;
            for (d=0;d<D;d++){
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[q-1+d*(L+1)] = j-1;
                B[q-1+d*(L+1)] = fabs(1.0-fabs(x[i+d*nx]-(j*2.0-1.0)/(mi-1))*(mi-1));
            }
        }

        for (ii=0;ii<nc;ii++){
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)];
            for (d=D-2;d>0;d--){
                b*=B[combs[ii+d*nc]-1+d*(L+1)];
                id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(L+1)];
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1])+covers_loc[ii]-1;
            for (d=0;d<nA;d++)
                y[i+d*nx] += b*w[id1+d*nG];
        }
    }
}



void interp_l(long D, long L, long nG, long nx, long nc, double* grid, double* x, double* w, double* Y, short* combs, short* combs_dM, int* covers_loc, int* mgrid)
{
    long i,ii,j,d,q,id1,dmi,mi;
    long J[D*(L+1)];
    double B[D*(L+1)],b,y;

    int pstart = 1;
    for (d =0;d<D;d++)
        pstart *=combs_dM[nc-1+d*nc];
    pstart+=covers_loc[nc-1]-1;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi,y)
    for (i=0;i<nx;i++){
        y=0.0;
        for (d=0;d<D;d++){
            J[d*(L+1)] = 0;
            J[d*(L+1)+1] = (x[i+d*nx]>0.5);
            B[d*(L+1)] = 1.0;
            B[d*(L+1)+1] = fabs(1.0-fabs(x[i+d*nx]-J[d*(L+1)+1])*(2));
        }
        for (q=3;q<=L+1;q++){
            dmi = pow(2,q-2);
            mi = dmi*2+1;
            for (d=0;d<D;d++){
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[q-1+d*(L+1)] = j-1;
                B[q-1+d*(L+1)] = fabs(1.0-fabs(x[i+d*nx]-(j*2.0-1.0)/(mi-1))*(mi-1));
            }
        }

        for (ii=0;ii<nc;ii++){
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)];
            for (d=D-2;d>0;d--){
                b*=B[combs[ii+d*nc]-1+d*(L+1)];
                id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(L+1)];
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1])+covers_loc[ii]-1;
            y+=b*w[id1];
        }

        for (ii=pstart;ii<nG;ii++){
            b=1.0;
            for (d=0;d<D;d++){
                b*=bf_cc_l(x[i+d*nx],grid[ii+d*nG],mgrid[ii+d*nG]);
                if (b==0.0)
                    break;
            }
            if (b!=0.0)
                y+=b*w[ii];
        }
        Y[i] = y;
    }
}


void interp_l(long D, long L, long nG, long nx, long nc, long nA, double* grid, double* x, double* w, double* Y, short* combs, short* combs_dM, int* covers_loc, int* mgrid)
{
    long i,ii,j,d,q,id1,dmi,mi;
    long J[D*(L+1)];
    double B[D*(L+1)],b;

    int pstart = 1;
    for (d =0;d<D;d++)
        pstart *=combs_dM[nc-1+d*nc];
    pstart+=covers_loc[nc-1]-1;

    #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi)
    for (i=0;i<nx;i++){
        for (d=0;d<D;d++){
            J[d*(L+1)] = 0;
            J[d*(L+1)+1] = (x[i+d*nx]>0.5);
            B[d*(L+1)] = 1.0;
            B[d*(L+1)+1] = fabs(1.0-fabs(x[i+d*nx]-J[d*(L+1)+1])*(2));
        }
        for (q=3;q<=L+1;q++){
            dmi = pow(2,q-2);
            mi = dmi*2+1;
            for (d=0;d<D;d++){
                j = round(x[i+d*nx]*dmi+0.5);
                if (j<1)
                    j=1;
                else if (j>dmi)
                    j = dmi;
                J[q-1+d*(L+1)] = j-1;
                B[q-1+d*(L+1)] = fabs(1.0-fabs(x[i+d*nx]-(j*2.0-1.0)/(mi-1))*(mi-1));
            }
        }

        for (ii=0;ii<nc;ii++){
            b = B[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)]*B[combs[ii]-1];
            id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(L+1)];
            for (d=D-2;d>0;d--){
                b*=B[combs[ii+d*nc]-1+d*(L+1)];
                id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(L+1)];
            }
            id1 = id1*combs_dM[ii]+(J[combs[ii]-1])+covers_loc[ii]-1;
            for (d=0;d<nA;d++)
                Y[i+d*nx] += b*w[id1+d*nG];
        }
        for (ii=pstart;ii<nG;ii++){
            b=1.0;
            for (d=0;d<D;d++){
                b*=bf_cc_l(x[i+d*nx],grid[ii+d*nG],mgrid[ii+d*nG]);
                if (b==0.0)
                    break;
            }
            if (b!=0.0){
                for (d=0;d<nA;d++)
                    Y[i+d*nx] += b*w[ii+d*nG];
            }
        }
    }
}
