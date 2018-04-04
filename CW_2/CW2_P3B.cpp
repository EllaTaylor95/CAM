#include <iostream>
#include <cmath>
#include <iomanip>
#include <math.h>

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

double normal_CFD(double x)
{
   return 0.5 * erfc(-x/sqrt(2));
}

double u_exact(double x, double t, double r, double sigma_sq, double K)
{
return K*exp(-r*t)*(normal_CFD(-pow(sigma_sq*t,-0.5)*(log(x/K)+(r-sigma_sq/2)*t)))-x*normal_CFD(-pow(sigma_sq*t,-0.5)*(log(x/K)+(r+sigma_sq/2)*t));
}

double g(double x, double K)
{
double g;
if(x<K)
    {
    g = K-x;
    }
else
    {
    g = 0;
    }
return g;
}

double f_0(double t, double K, double r)
{
return K*exp(-r*t);
}

//Assign vectors according to discretised Black Scholes Eq/
void assign_vectors(int m, int n, double K, double R, double r, double sigma_sq, double delta_t, double* x, double* f, double* U, double* diag, double* upper, double* lower, double* U_prev)
    {
double h = R/(double)(m);

f[0] = f_0((double)(n)*delta_t,K,r)*(-delta_t*sigma_sq*pow(0,2)/(2*pow(h,2)));
f[m-1] = u_exact(R,(double)(n)*delta_t,r,sigma_sq,K)*(-delta_t*r*R/h - delta_t*sigma_sq*pow(R,2)/(2*pow(h,2)));

for (int i=0; i<m; i++)
    {
    diag[i] = 1 + delta_t*r + delta_t*r*x[i]/h + delta_t*sigma_sq*pow(x[i],2)/pow(h,2);
    U_prev[i] = U[i] + f[i];
    }

for (int i=0; i<m-1; i++)
    {
    lower[i] = - delta_t*sigma_sq*pow(x[i],2)/(2*pow(h,2));
    upper[i] = - delta_t*r*x[i]/h - delta_t*sigma_sq*pow(x[i],2)/(2*pow(h,2));
    }

    }

int main()
{
//m x m is size of matrix A
int m = 10, num_of_time_steps = 5;
double R = 300.0, h = R/(double)(m), T = 5.0, delta_t = T/(double)(num_of_time_steps), r = 0.0, sigma_sq = 0.5, K = 100;
//Lower diagonal
double* lower;
lower = new double[m-1];
//Leading diagonal
double* diag;
diag = new double[m];
//Upper diagonal
double* upper;
upper = new double[m-1];
//Initial or previous value of U
double* U_prev;
U_prev = new double[m];
//Value of U being calculated
double* U;
U = new double[m];
//Boundary condition when t = 0 (initial conditions)
double* f;
f = new double[m];

double* x;
x = new double[m];


for(int i =0; i<m; i++)
    {
    //Assigning ICs
    x[i] = (i+1)*h; //Shifted across by one as U(0) already known
    f[i] = 0;
    U[i] = g(x[i],K);
    }

assign_vectors(m,0,K,R,r,sigma_sq,delta_t,x,f,U,diag,upper,lower,U_prev);

for(int n = 0; n<num_of_time_steps+1; n++)
    {
    tridiagonal_matrix_solver(m,U,lower,diag,upper,U_prev);

    std::cout<<std::endl<<std::endl<<std::setw(1)<<"i"<<std::setw(15)<<"U"<<std::endl;
    std::cout<<std::setw(1)<<0<<std::setw(15)<<f_0((double)(n)*delta_t,K,r)<<std::endl;
    for(int i = 0; i<m; i++)
        {
        std::cout<<std::setw(1)<<i+1<<std::setw(15)<<U[i]<<std::endl;
        }
    std::cout<<std::setw(1)<<m+1<<std::setw(15)<<u_exact(R,(double)(n)*delta_t,r,sigma_sq,K)<<std::endl;

    assign_vectors(m,n,K,R,r,sigma_sq,delta_t,x,f,U,diag,upper,lower,U_prev);
    }

//Define error as max (e_N) = max error at final time

double Ei, E_max = 0;

for(int i =0; i<m; i++)
    {
    Ei = fabs(U[i] - u_exact(x[i],T,r,sigma_sq,K));

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


