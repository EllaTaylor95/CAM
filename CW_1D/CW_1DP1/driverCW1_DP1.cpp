#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>

void funct(double* f, double* y, double alpha, double beta)
{
f[0] = y[1];
f[1] = -pow(alpha,2)*sin(y[0])-pow(beta,2)*y[1];
}

void Euler(int dim, double* yn, double* y_initial, double h, int n_max, double alpha, double beta)
{
for (int j = 0; j<dim; j++)
    {
    yn[j] = y_initial[j];
    }

double f[dim];
int i = 0;
double Energy;

std::cout<<std::setw(5)<<"tn"<<std::setw(15)<<"y"<<std::setw(25)<<"y'"<<std::setw(30)<<"Energy"<<std::endl;

while(i<=n_max)
    {

    Energy = 0.5*pow(yn[1],2)+pow(alpha,2)*(1-cos(yn[0]));

    std::cout<<std::setw(5)<<i*h<<std::setw(15)<<yn[0]<<std::setw(25)<<yn[1]<<std::setw(30)<<Energy<<std::endl;

    funct(f,yn,alpha, beta);

    for (int j = 0; j<dim; j++)
        {
        yn[j] = yn[j] + h*f[j];
        }
    i++;
    }
}

int main()
{
const double alpha = 2.0;
const double beta = 1.0;
const int dim = 2;
const double pi = 3.14159265359;

double* y;
y = new double[dim];
y[0] = pi/2.0;
y[1] = 0.0;

int n_max = 32;
double h = 0.25;

Euler(dim, y,y, h, n_max, alpha,beta);

return 0;
}
