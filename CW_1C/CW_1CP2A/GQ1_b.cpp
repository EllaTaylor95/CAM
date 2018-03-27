#include <cmath>
#include <iostream>

double f2(double x);

double CGQ1(int n, double a, double b)
{
double h = (b-a)/(double)(n);
double sum = 0;
double Xi;

//f = (tanh(100(x-2)))((3-x)/10)sqrt(x^2 - x sqrt(2)+ 0.5005)

for(double i = a; i<b; i = i+h)
    {
    Xi = (2.0 * i + h)/2.0;
    sum = sum + h * f2(Xi);
    }

return sum;
}
