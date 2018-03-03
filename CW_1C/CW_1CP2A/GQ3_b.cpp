#include <cmath>
#include <cassert>

double f2(double x);

//n i number of intervals (not number of  nodes)
double CGQ3(int n, double a, double b)
{
double sum = 0.0;
double h = (b-a)/((double)(n));

double x[3];
x[0] = -sqrt(3.0/5.0);
x[1] = 0;
x[2] = sqrt(3.0/5.0);

double Xi[3];

//f = (tanh(100(x-2)))((3-x)/10)sqrt(x^2 - x sqrt(2)+ 0.5005)

for (double i = a; i<b; i = i+h)
    {

    for (int j=0; j<3; j++)
    {
    //mapping points on to interval [-1,1]
    Xi[j] = x[j]*h/2.0 + (2.0*i+h)/2.0;
    }

    sum = sum + (h/2.0)*( (5.0/9.0)*f2(Xi[0]) + (8.0/9.0)*f2(Xi[1]) + (5.0/9.0)*f2(Xi[2]));
    }

return sum;
}
