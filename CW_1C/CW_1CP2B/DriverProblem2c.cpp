#include <iostream>
#include <cmath>

double adapt_int(double a, double b, double tau);

int main()
{
const double a = 0.0;
const double b = 3.0;
const double tau = pow(10,-3);

std::cout<<"Gauss Quad 1 = " << adapt_int(a,b,tau)<<"\n";

return 0;
}
