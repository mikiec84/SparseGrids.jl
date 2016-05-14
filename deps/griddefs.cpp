inline double bf_cc_l (double x, double xij, int mi)
{
    double dx;
    if (mi==1)
        return 1.0;
    dx = 1.0-((mi-1.0)*fabs(x-xij));
    if (dx>0)
        return dx;
    else return 0.0;
}

inline double bf_cc_q (double x, double xij, int mi)
{
    double dx;
    if (mi==1)
        return 1.0;
    dx = 1.0-pow((mi-1.0)*(x-xij),2);
    if (dx>0)
        return dx;
    else return 0.0;
}


inline double bf_m_l (double x, double xij, int mi)
{
    if (fabs(x-xij)<1/(mi-1.0))
        return 1.-(mi-1.)*fabs(x-xij);
    else
        return 0.0;
}

inline double bf_m_q (double x, double xij, int mi)
{
    double dx;
    dx = 1.0-pow((mi-1.0)*(x-xij),2);
    if (dx>0)
        return dx;
    else
        return 0.0;
}

inline int M_cc(long i)
{
    if (i==1)
        return 1;
    else
        return pow(2,i-1)+1;
}

inline int dM_cc(long i)
{
    if (i<3)
        return pow(2,i-1);
    else
        return pow(2,i-2);
}

inline double dxi_cc(long i, long j)
{
    if (i==1)
        return 0.5;
    else if (i==2)
        if (j==1)
            return 0.0;
        else
            return 1.0;
    else
        return (j*2.0-1.0)/(M_cc(i)-1);
}

inline int indexshift_cc(long i, long j, long q)
{
    if (i==1)
        return (M_cc(q)-1)/2+1;
    else if (i==2)
        if (j==1)
            return 1;
        else
            return M_cc(q);
    else
        return ((M_cc(q)-1)*(j*2-1))/(M_cc(i)-1)+1;

}
