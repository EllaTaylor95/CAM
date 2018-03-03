#include <iostream>
#include <cmath>

double GQ1(double a, double b);
double GQ3(double a, double b);

int main()
{
const double a = 0.0;
const double b = 3.0;

double GQ_1;

GQ_1 = GQ1(a,b);

std::cout<<"GQ1 = "<<GQ_1<<"\n";

double GQ_3;

GQ_3 = GQ3(a,b);

std::cout<<"GQ3 = "<<GQ_3<<"\n";

return 0;
}
