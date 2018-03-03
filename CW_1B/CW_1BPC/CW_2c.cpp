#include <iostream>
#include <cmath>
#include <cassert>

double evaluate_qh(int n, double h, double x_val, double x0, double* c)
{
double qh = 0;

//For any x in [xk-1,xk] bi is only locally non zero
//therefore q(x) = ck-2Bk-2(x) + ck-1Bk-1(x) + ckBk(x) + ck+1Bk+1(x)
int k = int((x_val - x0)/h +1.0);

assert(k>=0);

if (k<2)
    {
    qh = 4.0*c[k]+ 2.0*c[k+1];
    }
else if ((k>1) && (k<=n))
    {
    qh = c[k-1] +4.0*c[k]+c[k+1];
    }
else
    {
    qh = 2.0*c[k-1]+4.0*c[k];
    }
return qh;
}
