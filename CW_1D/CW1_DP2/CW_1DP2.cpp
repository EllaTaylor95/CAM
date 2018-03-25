#include <iostream>
#include <cmath>
#include <iomanip>

void funct(double* f, double* y, double alpha, double beta)
{
f[0] = y[1];
f[1] = -pow(alpha,2)*sin(y[0])-pow(beta,2)*y[1];
}

void Newton_Method(double TOL, double* y_1, double* y_2, int n_max, double alpha, double beta, double h)
{
double* u;
u = new double[2];

int counter = 0;

double *f;
f = new double[2];
double F[2];

double F_norm = 10.0;

for (int i = 1; i < n_max; i++)
    {
    //Setting initial guess y_n+1(0) = y_n
    y_1[i] = y_1[i-1];
    y_2[i] = y_2[i-1];

    while (F_norm >= TOL)
        {
        F_norm = 0;

        std::cout<<"counter = "<<counter<<"\n";
        counter++;

        u[0] = y_1[i];
        u[1] = y_2[i];

        funct(f,u,alpha,beta);

        F[0] = y_1[i] - y_1[i-1] - h*f[0];
        //std::cout<<"F[0] = "<<F[0]<<"\n";
        F[1] = y_2[i] - y_2[i-1] - h*f[1];
        //std::cout<<"F[1] = "<<F[1]<<"\n";

        y_1[i] = y_1[i] - ( 1.0/(1.0 + h*pow(alpha,2) * cos(y_1[i])) ) * (F[0] + h*F[1]);
       // std::cout<<"y_1["<<i<<"] = "<<y_1[i]<<"\n";
        y_2[i] = y_2[i] - ( 1.0/(1.0 + h*pow(alpha,2) * cos(y_1[i])) ) * (-h*pow(alpha,2)*cos(y_1[i])*F[0]+F[1]);
        //std::cout<<"y_2["<<i<<"] = "<<y_2[i]<<"\n";

        F_norm = pow(F[0],2) + pow(F[1],2);
        std::cout<<"F_norm = "<<F_norm<<"\n";
        }
    }
delete[] f;
delete[] u;
}

int main()
{
double alpha = 2.0;
double beta = 0.0;
double pi = 3.14159265359;

double T = 8.0;
int n_max = 32;
double h = T/(double)(n_max);

double* y_1;
double* y_2;

y_1 = new double[n_max];
y_2 = new double[n_max];


y_1[0] = pi/2.0;
y_2[0] = 0.0;

for (int i = 1; i<n_max; i++)
{
    y_1[i] = 0;
    y_2[i] = 0;
}

double TOL = pow(10,-12);

//From backward Euler have F(y_n+1) = y_n+1 - y_n -hf(y_n+1,t+1) = 0
//Then apply Newtons method to this

Newton_Method(TOL, y_1, y_2, n_max, alpha, beta, h);

delete[] y_1;
delete[] y_2;

return 0;
}
