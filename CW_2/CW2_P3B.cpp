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
   return 0.5 * erfc(-x/sqrt(2.0));
}

double d(double t, double x, double sign)
{
    double K =100.0, r = 0.0, sigma_sq = pow(0.5,2);

    return pow(sigma_sq*t,-0.5)*(log(x/K) + (r + sign*sigma_sq/2)*t);
}

double u_exact(double x, double t)
{
double K =100.0, r = 0.0, sigma_sq = pow(0.5,2);

return K*exp(-r*t)*(normal_CFD(-d(t,x,-1)))-x*normal_CFD(-d(t,x,1));
}

double g(double x)
{
double K = 100;
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

double f_0(double t)
{
double K = 100.0, r = 0.0;
return K*exp(-r*t);
}

double f_R(double t)
{
double R = 300.0;
return u_exact(R,t);
}

double alpha(int i, double delta_t)
{
double r = 0.0, sigma_sq = pow(0.5,2);
return  1.0 + delta_t*r + delta_t*r*i + delta_t*sigma_sq*pow(i,2);
}

double beta(int i, double delta_t)
{
double r = 0.0, sigma_sq = pow(0.5,2);
return -delta_t*r*i - delta_t*sigma_sq*pow(i,2)/2.0;
}

double gamma(int i, double delta_t)
{
double r = 0.0, sigma_sq = pow(0.5,2);
return -delta_t*sigma_sq*pow(i,2)/2.0;
}


//Assign vectors according to discretised Black Scholes Eq/
void assign_vectors(int m, double h, int n,double delta_t, double* f, double* U, double* diag, double* upper, double* lower, double* U_prev)
{
f[0] = f_0((double)(n)*delta_t)*gamma(1.0,delta_t);
f[m-1] = f_R((double)(n)*delta_t)*beta((double)(m),delta_t);

for (int i=0; i<m; i++)
    {
    diag[i] = alpha((double)(i+1),delta_t);
    U_prev[i] = U[i] - f[i];
    }

for (int i=0; i<m-1; i++)
    {
    lower[i] = gamma((double)(i+2),delta_t);
    upper[i] = beta((double)(i+1),delta_t);
    }

    }


double E_max(int m, double h, double t, double*U)
{
//Define error as max (e_N) = max error at final time
double Ei, E_max = 0;

for(int i =0; i<m; i++)
    {
    Ei = fabs(U[i] - u_exact(h*(i+1),t));

    if(E_max < Ei)
        {
        E_max = Ei;
        }
    }
return E_max;
}

int main()
{
//m x m is size of matrix A
int m = 19, num_of_time_steps = 10; //m is number of interior nodes
double R = 300.0, h = R/((double)(m)+1), T = 5.0, delta_t = T/(double)(num_of_time_steps);

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
double* f;
f = new double[m];

for(int i =0; i<m; i++)
    {
    //Assigning ICs
    f[i] = 0;
    U[i] = g(h*(double)(i+1));
    }

std::cout<<"n = 0 \n";
std::cout<<std::setw(1)<<"i"<<std::setw(15)<<"U"<<std::setw(25)<<"u(x_i)"<<std::endl;
std::cout<<std::setw(1)<<0<<std::setw(15)<<f_0((double)(0)*delta_t)<<std::setw(25)<<u_exact(0,(double)(0)*delta_t)<<std::endl;
for(int i = 0; i<m; i++)
        {
    std::cout<<std::setw(1)<<i+1<<std::setw(15)<<U[i]<<std::setw(25)<<u_exact((i+1)*h,(double)(0)*delta_t)<<std::endl;
        }
std::cout<<std::setw(1)<<m+1<<std::setw(15)<<f_R((double)(0)*delta_t)<<std::setw(25)<<u_exact(R,(double)(0)*delta_t)<<std::endl;
std::cout<<"E_max = "<<E_max(m,h,0*delta_t,U)<<std::endl;

for(int n = 1; n<=num_of_time_steps; n++)
    {
    assign_vectors(m,h,n,delta_t,f,U,diag,upper,lower,U_prev);
    tridiagonal_matrix_solver(m,U,lower,diag,upper,U_prev);
    std::cout<<"\n"<<"n = "<<n<<std::endl;
    std::cout<<std::setw(1)<<"i"<<std::setw(15)<<"U"<<std::setw(25)<<"u(x_i)"<<std::endl;
    std::cout<<std::setw(1)<<0<<std::setw(15)<<f_0((double)(n)*delta_t)<<std::setw(25)<<u_exact(0,(double)(n)*delta_t)<<std::endl;
    for(int i = 0; i<m; i++)
        {
        std::cout<<std::setw(1)<<i+1<<std::setw(15)<<U[i]<<std::setw(25)<<u_exact((i+1)*h,(double)(n)*delta_t)<<std::endl;
        }
    std::cout<<std::setw(1)<<m+1<<std::setw(15)<<f_R((double)(n)*delta_t)<<std::setw(25)<<u_exact(R,(double)(n)*delta_t)<<std::endl;
    std::cout<<"E_max = "<<E_max(m,h,(double)(n)*delta_t,U)<<std::endl;

    }

delete[] lower;
delete[] upper;
delete[] diag;
delete[] U_prev;
delete[] U;

return 0;
}


