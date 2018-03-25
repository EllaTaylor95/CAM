#include <iostream>
#include <cmath>

void tridiagonal_matrix_solver(int n, double* U, double* lower, double* diag, double* upper, double* f)
    {
    //Elimination stage

    //f[0] = f[0] and d[0] = d[0] no change
    for (int i = 1; i<n; i++)
        {
        diag[i] = diag[i] - ((upper[i-1]*lower[i-1])/diag[i-1]);
        f[i]    = f[i] - ((f[i-1]*lower[i-1])/diag[i-1]);
        }

    //Backsolving
    //Bottom row is a special case
    U[n-1] = f[n-1]/diag[n-1];

    for (int i = n-2; i >= 0; i--)
        {
        U[i] = (f[i] - upper[i]*U[i+1])/diag[i];
        }

    }

int main()
{
//m x m is size of matrix A
int m = 20;
double h = 1/(double)(m);
double pi = 3.14159265359;

double* lower;
lower = new double[m-1];

double* diag;
diag = new double[m];

double* upper;
upper = new double[m-1];

double* f;
f = new double[m];


double* U;
U = new double[m];


U[0] = 0.0;
diag[0] = 2.0/pow(h,2);
upper[0] = -1.0/pow(h,2);
//f(x0) + alpha
f[0] = 0.0 ;
lower[m-2] = -1.0/pow(h,2);
diag[m-1] = 2.0/pow(h,2);
//f(xm) + beta
f[m-1] = -pow(pi,2)*sin(1.0*pi) + 1.0/pow(h,2);
U[m-1] = 0;

for (int i=1; i<m-1; i++)
    {
    U[i] = 0;
    lower[i-1] = -1.0/pow(h,2);
    diag[i] = 2/pow(h,2);
    upper[i] = -1/pow(h,2);
    f[i] = -pow(pi,2)*sin(pi*i*h);
    }

tridiagonal_matrix_solver(m,U,lower,diag,upper,f);

double U_last = U[m-1];

for(int i = m-1; i>0; i--)
    {
    U[i] = U[i-1];
    }
U[0] = 0;

for(int i = 0; i<m; i++)
    {
    std::cout<<"U["<<i<<"] = "<<U[i]<<"\n";
    }
std::cout<<"U["<<m<<"] = "<<U_last<<"\n"<<"U["<<m+1<<"] = "<<1<<"\n";


delete[] lower;
delete[] upper;
delete[] diag;
delete[] f;
delete[] U;

return 0;
}

