#include <iostream>
#include <cmath>
#include <iomanip>

double alpha(int i)
{
int m = 16, num_of_iterations = 20;
double R = 300.0, h = R/(double)(m), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25;

double alpha = 1.0 + r*delta_t*(1.0+i) + sigma_sq*pow(i,2)*delta_t;

return alpha;
}

double beta(int i)
{
int m = 16,num_of_iterations = 20;
double R = 300.0, h = R/(double)(m), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25;

double beta = -0.5*delta_t*sigma_sq*pow(i,2) - r*(i)*delta_t;

return beta;
}

double gamma(int i)
{
int m = 16,num_of_iterations = 20;
double R = 300.0, h = R/(double)(m), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25;

double gamma = -0.5*delta_t*sigma_sq*pow(i,2);

return gamma;
}

int main()
{

int m = 16,num_of_iterations = 20;
double R = 300.0, h = R/(double)(m), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25, K = 100;

double u_bar[m+1], u[m+1], g[m+1];

for (int i = 0; i<m+1; i++)
    {
    u[i] = 0.0;

    if(i*h<K)
        {
        g[i] = K - i*h;
        }
    else
        {
        g[i] = 0.0;
        }
    }

u[0] = K;
u[m] = 0.0;

int counter =0;

while(counter<num_of_iterations)
    {

    std::cout<<"n"<<std::setw(4)<<"u"<<std::endl;

    for(int i = 1; i<m; i++)
        {
        u_bar[i] = (1.0/alpha(i))*( u[i] - gamma(i-1)*u[i-1] - beta(i+1)*u[i+1] );
        }

    for(int i = 1; i<m; i++)
        {

        if (g[i] >= (omega*u_bar[i]+(1.0-omega)*u[i]) )
            {
            u[i] = g[i];
            }
        else
            {
            u[i] = omega*u_bar[i]+(1-omega)*u[i];
            }
        }
    counter ++;

    std::cout<<counter<<std::endl;
    for (int i = 0; i<m+1; i++)
        {
        std::cout<<std::setw(5)<<u[i]<<std::endl;
        }
    }

return 0;
}
