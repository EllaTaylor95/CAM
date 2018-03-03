#include <cmath>


double f1(double x)
{
const double pi = 3.14159265358979323846;
return (1.0/18.0)*pow(x,4) - sin(pi*x/6.0);
}


double f2(double x)
{
return tanh(100.0*(x-2.0)) * ((3.0-x)/10.0) * pow(pow(x,2)-x*pow(2,0.5)+0.5005,0.5);
}
