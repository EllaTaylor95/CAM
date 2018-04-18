#include <iostream>
#include <cmath>
#include <iomanip>

int main()
{
//mesh size
int m = 16;
double h = 1.0/(double)(m), phi = 1.0, f = 50/3, omega = 1.8;

double u_bar[m], u[m], u_prev[m];

for (int i = 0; i<m; i++)
    {
    u[i] = 0.0;
    }

int counter =0;

while(counter<8)
    {

    std::cout<<"n"<<std::setw(4)<<"u"<<std::endl;

    u_bar[0] = u[1]*0.5;

    for(int i = 1; i<m; i++)
        {
        u_bar[i] = 0.5*( f*pow(h,2) + u[i-1] + u_prev[i+1] );
        }

    for(int i = 0; i<m; i++)
        {

        u_prev[i] = u[i];

        if (phi <= (omega*u_bar[i]+(1-omega)*u_prev[i]) )
            {
            u[i] = phi;
            }
        else
            {
            u[i] = omega*u_bar[i]+(1-omega)*u_prev[i];
            }
        }
    counter ++;

    std::cout<<counter<<std::endl;
    for (int i = 0; i<m; i++)
        {
        std::cout<<std::setw(5)<<u[i]<<std::endl;
        }
    }



return 0;
}
