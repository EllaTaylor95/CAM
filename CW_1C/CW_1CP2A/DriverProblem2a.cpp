#include <iostream>
#include <iomanip>
#include <cmath>

double CGQ1(int n, double a, double b);
double CGQ3(int n, double a, double b);

int main()
{
int n; //n is number of intervals
double h;
const double a = 0;
const double b = 3;

std::cout<< "N"<<std::setw(20)<<"h"<<std::setw(30)<<"CGQ1"<<std::setw(35)<<"CGQ3"<<std::endl;
std::cout<<"------------------------------------------------------------------------------------------"<<std::endl;

for (int i=0; i<=7; i++)
{
n = pow(2,i);
h = (b-a)/n;
std::cout<<n<<std::setw(20)<<h<<std::setw(30)<<CGQ1(n,a,b)<<std::setw(35)<<CGQ3(n,a,b)<<std::endl;
}

return 0;
}
