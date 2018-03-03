#include <iostream>
#include <ctime>

void tridiagonal_matrix_solver(int n, double*x, double* lower, double* diag, double* upper, double* f);

int main()
{
int n;
std::cout<< "Enter an integer n \n";
std::cin>>n;

double* lower;
lower = new double[n-1];

double* diag;
diag = new double[n];

double* upper;
upper = new double[n-1];

double* f;
f = new double[n];

double* x;
x = new double[n];

x[0] = 0;
diag[0] = 4;
upper[0] = 2;
f[0] = 1;
lower[n-2] = 2;
diag[n-1] = 4;
f[n-1] = (double)(n);
x[n-1] = 0;

for (int i=1; i<n-1; i++)
    {
    x[i] = 0;
    lower[i-1] = 1;
    diag[i] = 4;
    upper[i] = 1;
    f[i] = (double)(i)+1;
    }

clock_t start_s = clock();
tridiagonal_matrix_solver(n,x,lower,diag,upper,f);
clock_t stop_s = clock();
std::cout<<"Time = "<<(stop_s - start_s)/((double)(CLOCKS_PER_SEC))<<"\n";

std::cout<<"x[1] = "<<x[0]<<"\n";
std::cout<<"x[n] = "<<x[n-1]<<"\n";

delete[] lower;
delete[] upper;
delete[] diag;
delete[] f;
delete[] x;

return 0;
}
