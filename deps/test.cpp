//  void cinterp(long D, long Q, long nx, long nc,
//     double* x, double* w, double* y,
//     short* combs, short* combs_dM, int* combl,
//     long* topc, uint32_t* Tp, long* Ts)
// {
//     long i,ii,iii,j,d,q,id1,dmi,mi;
//     long J[D*(Q+1)],id[D], cTs[D];
//     size_t nT;
//     double B[D*(Q+1)],b;
//     uint32_t *t0,*t1;
//     char Todd;
//     t0 = new uint32_t[Ts[0]];
//     t1 = new uint32_t[Ts[0]];
//
//     for (i=0;i<Ts[0];i++){
//       t0[i] = 1;
//       t1[i] = 1;
//     }
//
//     i =0;
//     for (d=0;d<D;d++){
//         i+=Ts[d]*Ts[D+d];
//         cTs[d] = i;
//     }
//
//     #pragma omp parallel for private(d,q,j,J,B,ii,b,id1,dmi,mi,t0,t1,Todd,id,nT)
//     for (i=0;i<nx;i++)
//     {
//         for (d=0;d<D;d++){
//             J[d*(Q+1)]   = 0;
//             J[d*(Q+1)+1] = (x[i+d*nx]>0.5);
//             B[d*(Q+1)]   = 1.0;
//             B[d*(Q+1)+1] = fabs(1.0-fabs(x[i+d*nx]-J[d*(Q+1)+1])*(2));
//         }
//         for (q=3;q<=Q+1;q++)
//         {
//             dmi = pow(2,q-2);
//             mi = dmi*2+1;
//             for (d=0;d<D;d++)
//             {
//                 j = round(x[i+d*nx]*dmi+0.5);
//                 if (j<1)
//                         j=1;
//                 else if (j>dmi)
//                     j = dmi;
//                 J[q-1+d*(Q+1)] = j-1;
//                 B[q-1+d*(Q+1)] = fabs(1.0-fabs(x[i+d*nx]-(j*2.0-1.0)/(mi-1))*(mi-1));
//             }
//         }
//
//
//         for (iii=0;iii<D;iii++)
//         {
//             ii = topc[iii]-1;
//             id1 = J[combs[ii+(D-1)*nc]-1+(D-1)*(Q+1)];
//             for (d=D-2;d>0;d--)
//             {
//                 id1 = id1*combs_dM[ii+d*nc] + J[combs[ii+(d)*nc]-1+(d)*(Q+1)];
//             }
//             id[iii] = id1*combs_dM[ii]+(J[combs[ii]-1]);
//         }
//
//         nT = intersect(Tp + id[0]*Ts[0],Ts[0],Tp + id[1]*Ts[1] + cTs[0],Ts[1],t0);
//         Todd = false;
//
//         // for (d=1;d<D-1;d++){
//         //     if (!Todd){
//         //         nT = intersect(Tp+(cTs[d]+id[d+1]*Ts[d+1]),Ts[d+1],t0,nT,t1);
//         //     }
//         //     else{
//         //         nT = intersect(Tp+(cTs[d]+id[d+1]*Ts[d+1]),Ts[d+1],t1,nT,t0);
//         //     }
//         //     Todd = !Todd;
//         // }
//
//
//         for (ii=0;ii<nc;ii++){
//             b = 1.0;
//             for (d=0;d<D;d++)
//                 b*=B[combs[ii+d*nc]-1+d*(Q+1)];
//             // y[i] += b*w[t0[ii]-1];
//         }
//     }
//     delete[] t0;
//     delete[] t1;
// }
//
//
//
// void interp_l_t(int D, int Q, int nG, int nx,
//                 double * grid, short * lvl_s,
//                 double * w,
//                 double * xi, double * y,
//                 int * nextid,
//                 int is, int ie)
// {
//     double temp2;
//     int i, ii, d;
//
//     for (i=is;i<ie;i++)
//     {
//         ii=1;
//         while (ii<=nG)
//         {
//             temp2 =1.0;
//             for (d=0;d<D;d++)
//             {
//                 temp2 *= bf_cc_l( xi[i+d*nx] , grid[ii-1+d*nG] , lvl_s[ii-1+d*nG] );
//                 if (temp2==0.0)
//                     break;
//             }
//             if (temp2==0.0)
//                 ii+=1;
//             else
//             {
//                 y[i] += temp2*w[ii-1];
//                 ii=nextid[ii-1];
//             }
//         }
//     }
// }
//
// void interp_l_main(int D, int Q, int nG, int nx,
//                 double * grid, short * lvl_s,
//                 double * w,
//                 double * xi, double * y,
//                 int * nextid){
//   int nthreads = 8;
//   int dnx = nx/8;
//   interp_l_t(D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,7*dnx,nx);
//   thread t1(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,0*dnx,1*dnx);
//   thread t2(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,1*dnx,2*dnx);
//   thread t3(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,2*dnx,3*dnx);
//   thread t4(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,3*dnx,4*dnx);
//   thread t5(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,4*dnx,5*dnx);
//   thread t6(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,5*dnx,6*dnx);
//   thread t7(interp_l_t,D,Q,nG,nx,grid,lvl_s,w,xi,y,nextid,6*dnx,7*dnx);
//   t1.join();
//   t2.join();
//   t3.join();
//   t4.join();
//   t5.join();
//   t6.join();
//   t7.join();
// }
