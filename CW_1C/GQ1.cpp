#include <cmath>
#include <cassert>

double f1(double x);

double GQ1(double a, double b)
{
assert(a < b);
//Mapping Xi on to x using x = (b-a)/2 * Xi + (a+b)/2
//Integral approx equal to 2*f(0) for integrating between -1,1.
// Xi = 0 the x = (a+b)/2
//GQ1 = 2(b-a)/2 *f((a+b)/2)= (b-a)*f((a+b)/2)

return (b-a)*f1((a+b)/2);
}
