#include <cmath>
#include <cassert>

double f1(double x);

double GQ3(double a, double b)
{
assert(a < b);
//Mapping Xi on to x using x = (b-a)/2 * Xi + (a+b)/2
//GQ3 = 5/9 * f(-sqrt(3/5)) + 8/9 * f(0) + 5/9  * f(sqrt(3/5))
//f(x) = 1/18 * pow(x,4) - sin(pi x /6)
//sin functions cancel as odd
double x[3];
x[0] = -sqrt(3.0/5.0);
x[1] = 0;
x[2] = sqrt(3.0/5.0);

double Xi[3];
for (int i=0; i<3; i++)
    {
    Xi[i] = x[i]*(b-a)/2.0 + (a+b)/2.0;
    }

return ((b-a)/2.0) * ( (5.0/9.0) * f1(Xi[0]) +(8.0/9.0) * f1(Xi[1]) + (5.0/9.0) * f1(Xi[2]) );

}
