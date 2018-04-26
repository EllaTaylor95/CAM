#include <iostream>
#include <cmath>
#include <iomanip>

int main()
{
//mesh size
int m = 15;
double h = 1.0/((double)(m)+1), phi = 1.0, f = 50/3, omega = 1.8;

double u_bar[m], u[m], u_prev[m];

for (int i = 0; i<m; i++)
    {
    u[i] = 0.0;
    }

int counter =0;
double E_max = 1.0, ei;

while(E_max > pow(10, -5))
    {
    E_max = 0.0;
    u_bar[0] = 0.5*(f*pow(h,2) + u[1]);

    if (phi <= (omega*u_bar[0]+(1-omega)*u[0]) )
        {
        u[0] = phi;
        }
    else
        {
        u[0] = omega*u_bar[0]+(1-omega)*u[0];
        }


    for(int i = 1; i<m-1; i++)
        {
        u_bar[i] = 0.5*( f*pow(h,2) + u[i-1] + u[i+1] );

        if (phi <= (omega*u_bar[i]+(1-omega)*u[i]) )
            {
            u[i] = phi;
            }
        else
            {
            u[i] = omega*u_bar[i]+(1-omega)*u[i];
            }
        }

    u_bar[m-1] = 0.5*(f*pow(h,2) + u[m-2]);

    if (phi <= (omega*u_bar[m-1]+(1-omega)*u[m-1]) )
        {
        u[m-1] = phi;
        }
    else
        {
        u[m-1] = omega*u_bar[m-1]+(1-omega)*u[m-1];
        }

    counter ++;

    for (int i = 0; i<m; i++)
        {
        ei = fabs(u[i] - u_prev[i]);
        if(ei>E_max)
            {
            E_max = ei;
            }
        u_prev[i] = u[i];
        }

    std::cout<<"Iteration "<<counter<<std::endl;
    std::cout<<std::setw(5)<<"i"<<std::setw(15)<<"U"<<std::endl;
    std::cout<<std::setw(5)<<0<<std::setw(15)<<0<<std::endl;

    for (int i = 0; i<m; i++)
        {
        std::cout<<std::setw(5)<<i+1<<std::setw(15)<<u[i]<<std::endl;
        }
    std::cout<<std::setw(5)<<m+1<<std::setw(15)<<0<<std::endl<<std::endl;
    }



return 0;
}
