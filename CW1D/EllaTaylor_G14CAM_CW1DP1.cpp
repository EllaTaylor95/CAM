#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>

void funct(double* f, double* y, double alpha, double beta)
{
f[0] = y[1];
f[1] = -pow(alpha,2)*sin(y[0])-pow(beta,2)*y[1];
}

void Forward_Euler(double* y, double h, int n_max, double alpha, double beta)
{

double f[2];
int i = 0;
double Energy;

std::cout<<std::setw(5)<<"tn"<<std::setw(15)<<"y1"<<std::setw(25)<<"y2"<<std::setw(30)<<"Energy"<<std::endl;

while(i<=n_max)
    {

    Energy = 0.5*pow(y[1],2)+pow(alpha,2)*(1-cos(y[0]));

    std::cout<<std::setw(5)<<i*h<<std::setw(15)<<y[0]<<std::setw(25)<<y[1]<<std::setw(30)<<Energy<<std::endl;

    funct(f,y,alpha, beta);

    for (int j = 0; j<2; j++)
        {
        y[j] = y[j] + h*f[j];
        }
    i++;
    }
}

int main()
{
const double alpha = 2.0;
const double beta = 0.1;
const double pi = 3.14159265359;

double* y;
y = new double[2];
y[0] = pi/2.0;
y[1] = 0.0;

int n_max = 100000;
//double T = 32;
double h = 0.0034;

Forward_Euler(y, h, n_max, alpha,beta);

return 0;
}
