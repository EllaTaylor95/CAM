#include <iostream>
#include <cmath>
#include <iomanip>

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


void assign_vectors(int m, double delta_t, double a, double* g, double* U, double* diag, double* upper, double* lower, double* U_prev)
    {
double h = 1/(double)(m);
double pi = 3.14159265359;

for (int i=0; i<m; i++)
    {
    diag[i] = 1.0+ 2.0*a*delta_t/pow(h,2);
    U_prev[i] = U[i] + g[i];
    }

for (int i=0; i<m-1; i++)
    {
    lower[i] = -a*delta_t/pow(h,2);
    upper[i] = -a*delta_t/pow(h,2);
    }


    }

int main()
{
//m x m is size of matrix A
int m = 10, num_time_steps = 5;
double T = 5, h = 1/(double)(m), pi = 3.14159265359, a = pow(pi,-2), delta_t = T/num_time_steps;

double* lower;
lower = new double[m-1];

double* diag;
diag = new double[m];

double* upper;
upper = new double[m-1];

double* U_prev;
U_prev = new double[m];

double* U;
U = new double[m];

double* g;
g = new double[m];

double h = 1/(double)(m);
double pi = 3.14159265359;
double a = pow(pi,-2);
double delta_t = 0.25;

for(int i =0; i<m; i++)
    {
    g[i] = 0;
    U[i] = sin(2.0*pi*h*(double)(i+1))+(double)(i+1)*h;
    }
g[m-1] = a*delta_t/pow(h,2);

assign_vectors(m,delta_t,a,g,U,diag,upper,lower,U_prev);

for(int n = 0; n<4; n++)
    {
    tridiagonal_matrix_solver(m,U,lower,diag,upper,U_prev);

    std::cout<<std::endl<<std::endl<<std::setw(1)<<"i"<<std::setw(15)<<"U"<<std::endl;
    std::cout<<std::setw(1)<<0<<std::setw(15)<<0<<std::endl;
    for(int i = 0; i<m; i++)
        {
        std::cout<<std::setw(1)<<i+1<<std::setw(15)<<U[i]<<std::endl;
        }
    std::cout<<std::setw(1)<<m+1<<std::setw(15)<<1<<std::endl;

    assign_vectors(m,delta_t,a,g,U,diag,upper,lower,U_prev);
    }

//Define error as max (e_N) = max error at final time

double Ei, E_max = 0;

for(int i =0; i<m; i++)
    {
    Ei = (double)(i+1)*h + exp(-4)*sin(pi*(double)(i+1)*h) - U[i];

    if(E_max < Ei)
        {
        E_max = Ei;
        }
    }

std::cout<<"E_max = "<<E_max<<std::endl;

delete[] lower;
delete[] upper;
delete[] diag;
delete[] U_prev;
delete[] U;

return 0;
}

