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


void reset_vectors(int m, double* U, double* diag, double* U_init)
    {
double h = 1/(double)(m);
double pi = 3.14159265359;
double a = 1.0;
double delta_t = 0.01;

for (int i=0; i<m; i++)
    {
    U[i] = sin(2.0*pi*h*(double)(i))+(double)(i)*h;
    diag[i] = 1.0+ 2.0*a*delta_t/pow(h,2);
    U_init[i] = U[i];
    U[i] = 0;
    }

    }

int main()
{
//m x m is size of matrix A
int m = 10;
double a = 1.0;
double delta_t = 0.01;
double h = 1/(double)(m);
double pi = 3.14159265359;

double* lower;
lower = new double[m-1];

double* diag;
diag = new double[m];

double* upper;
upper = new double[m-1];

double* U_init;
U_init = new double[m];

double* U;
U = new double[m];

for (int i=0; i<m-2; i++)
    {
    lower[i] = -a*delta_t/pow(h,2);
    upper[i] = -a*delta_t/pow(h,2);
    }

U_init[m-1] = U_init[m-1] + a*delta_t/pow(h,2);


for(int n = 0; n<4; n++)
    {
    n++;
    reset_vectors(m,U,diag,U_init);
    tridiagonal_matrix_solver(m,U,lower,diag,upper,U_init);

double U_last = U[m-1];

for(int i = m-1; i>0; i--)
    {
    U[i] = U[i-1];
    }
U[0] = 0;

std::cout<<std::endl<<std::endl<<std::setw(1)<<"i"<<std::setw(15)<<"U"<<std::endl;
for(int i = 0; i<m; i++)
    {
    std::cout<<std::setw(1)<<i<<std::setw(15)<<U[i]<<std::endl;
    }
std::cout<<std::setw(1)<<m<<std::setw(15)<<U_last<<std::endl<<std::setw(1)<<m+1<<std::setw(15)<<1<<std::endl<<std::endl;

    }



//Define error as ||E||2 = sqrt(h*sum(|Ej^2|)) Ej = u(xj) - Uj

//double Ei, E_2norm = 0;
//
//for(int i =0; i<m; i++)
//    {
//    Ei = i*h - sin(pi*i*h) - U[i];
//    E_2norm = E_2norm + pow(Ei,2);
//    }
//E_2norm = sqrt(h*E_2norm);
//
//std::cout<<m<<std::setw(7)<<h<<std::setw(15)<<E_2norm<<std::endl;

delete[] lower;
delete[] upper;
delete[] diag;
delete[] U_init;
delete[] U;
    //}
return 0;
}

