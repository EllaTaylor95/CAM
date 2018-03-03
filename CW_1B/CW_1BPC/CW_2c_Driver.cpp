#include <iostream>
#include <cmath>
#include <iomanip>

void tridiagonal_matrix_solver(int n, double* x, double* lower, double* diag, double* upper, double* f);
double evaluate_qh(int n, double h, double x_val, double x0, double* c);

int main()
{
double h; //stepsize
double L = 2; //Length of  interval
int n; //n is the number of nodes
double pi = 3.14159265358979323846;

double x_val;
std::cout<<"Enter x \n";
std::cin>>x_val;
double q_h;
double Err = 0;

std::cout << "--------------------------------------------------------------------------"
            << std::endl;

std::cout << std::setw(10) << "h"
              << std::setw(20) << "q_h"
              <<std::setw(40) << "Err"
              << std::endl
              << "--------------------------------------------------------------------------"
            << std::endl;
double* x;
double* lower;
double* diag;
double* upper;
double* f;
double* c;
double* c_end;

double f_diff[2];
f_diff[0] = exp(x[0])*(sin(5.0*pi*x[0]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[0]/4.0));

for (int j = 1; j<15; j++)
    {

    h = pow(0.5,j);
    n = int((L/h) +1.0);
    x = new double[n];
    lower = new double[n-1];
    diag = new double[n];
    upper = new double[n-1];
    f = new double[n];
    c = new double[n];
    c_end = new double[n+2]; //c including end values

    for (int i=0; i<n; i++)
        {
        x[i] = 0.0 + h*(double)(i); //Interval (0,l)
        f[i] = exp(x[i])*sin((5.0/4.0)*pi*x[i]); //f(x) = exp(x)*sin(5pi/4)x)
        }
    f[0]   = f[0] + (h/3.0)*exp(x[0])*(sin(5.0*pi*x[0]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[0]/4.0));
    f[n-1] = f[n-1] - (h/3.0)*exp(x[n-1])*(sin(5.0*pi*x[n-1]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[n-1]/4.0));

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

    f_diff[1] = exp(x[n-1])*(sin(5.0*pi*x[n-1]/4.0)+(5.0*pi/4.0)*cos(5.0*pi*x[n-1]/4.0));

    c_end[0] = c[2] - (1.0/3.0)*h*f_diff[0];
    c_end[n+1] = c[n-2]+(1.0/3.0)*h*f_diff[1];


    for( int i = 0; i<n; i++)
    {
    c_end[i+1] = c[i];
    }

    q_h = evaluate_qh(n,h,x_val,x[0],c_end);

    Err = fabs(exp(1.0/3.0)*sin(pi*5.0/12.0) - q_h);


    std::cout<< std::setw(10)<<h
    <<std::setw(20)<<q_h
    <<std::setw(40)<<Err<<std::endl;

    }

delete[] lower;
delete[] upper;
delete[] diag;
delete[] f;
delete[] c;
delete[] x;
delete[] c_end;
}


