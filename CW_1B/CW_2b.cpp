#include <iostream>
#include <cmath>
#include <iomanip>

void tridiagonal_matrix_solver(int n, double* c, double* lower, double* diag, double* upper, double* f);

int main()
{
double h; //stepsize
std::cout<<"Enter stepsize h \n";
std::cin>>h;
double L; //Length of  interval
std::cout<<"Enter length of interval (0,L) \n";
std::cin>>L;
int n; //n is the number of nodes
n = int((L/h) +1.0);

double x[n];
double pi = 3.14159265358979323846;

double* lower;
lower = new double[n-1];

double* diag;
diag = new double[n];

double* upper;
upper = new double[n-1];

double* f;
f = new double[n];

for (int i=0; i<n; i++)
    {
    x[i] = 0.0 + h*(double)(i); //Interval (0,l)
    f[i] = exp(x[i])*sin((5.0/4.0)*pi*x[i]); //f(x) = exp(x)*sin(5pi/4)x)
    }


f[0]   = f[0] + (h/3.0)*exp(x[0])*(sin(5.0*pi*x[0]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[0]/4.0));
f[n-1] = f[n-1] - (h/3.0)*exp(x[n-1])*(sin(5.0*pi*x[n-1]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[n-1]/4.0));

//f dash evaluated at end points
double f_diff[2];
f_diff[0] = exp(x[0])*(sin(5.0*pi*x[0]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[0]/4.0));
f_diff[1] = exp(x[n-1])*(sin(5.0*pi*x[n-1]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[n-1]/4.0));

double* c;
c = new double[n]; //c not including c[-1] and c[n+1]

c[0] = 0.0;
diag[0] = 4.0;
upper[0] = 2.0;
lower[n-2] = 2.0;
diag[n-1] = 4.0;
c[n-1] = 0.0;

for (int i=1; i<n-1; i++)
    {
    c[i] = 0.0;
    lower[i-1] = 1.0;
    diag[i] = 4.0;
    upper[i] = 1.0;
    }

tridiagonal_matrix_solver(n,c,lower,diag,upper,f);

//Printing out entries of c
double c_end[n+2]; //c including end values

c_end[0] = c[2] - (1.0/3.0)*h*f_diff[0];
c_end[n+1] = c[n-2]+(1.0/3.0)*h*f_diff[1];

std::cout<<"c[-1] = "<<c[0]<<"\n";
for( int i = 0; i<n; i++)
{
c_end[i+1] = c[i];
std::cout<<"c["<<i<<"] = "<<c_end[i+1]<<"\n";
}
std::cout<<"c["<<n<<"] = "<<c_end[n+1]<<"\n";

delete[] lower;
delete[] upper;
delete[] diag;
delete[] f;
delete[] c;

return 0;
}
