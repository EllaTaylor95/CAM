#include <iostream>
#include <cmath>
#include <iomanip>

double normal_CFD(double x)
{
   return 0.5 * erfc(-x/sqrt(2));
}

double u_exact(double x, double t)
{
double K =100.0, r = 0.0, sigma_sq = pow(0.5,2);

return K*exp(-r*t)*(normal_CFD(-pow(sigma_sq*t,-0.5)*(log(x/K)+(r-sigma_sq/2)*t)))-x*normal_CFD(-pow(sigma_sq*t,-0.5)*(log(x/K)+(r+sigma_sq/2)*t));
}

double alpha(int i)
{
int m = 15, num_of_iterations = 20;
double R = 300.0, h = R/((double)(m)+1), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25;

double alpha = 1.0 + r*delta_t*(1.0+i) + sigma_sq*pow(i,2)*delta_t;

return alpha;
}

double beta(int i)
{
int m = 15,num_of_iterations = 20;
double R = 300.0, h = R/((double)(m)+1), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25;

double beta = -0.5*delta_t*sigma_sq*pow(i,2) - r*(i)*delta_t;

return beta;
}

double gamma(int i)
{
int m = 15,num_of_iterations = 20;
double R = 300.0, h = R/((double)(m)+1), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25;

double gamma = -0.5*delta_t*sigma_sq*pow(i,2);

return gamma;
}

int main()
{

int m = 15,num_of_iterations = 20;
double R = 300.0, h = R/((double)(m)+1), T = 5.0, delta_t = T/(double)(num_of_iterations), omega = 1.8, r = 0.05, sigma_sq = 0.25, K = 100;

double u_bar[m], u[m], g[m];

for (int i = 0; i<m; i++)
    {

    if((i+1)*h<K)
        {
        g[i] = K - (i+1)*h;
        }
    else
        {
        g[i] = 0.0;
        }
    u[i] = g[i];
    }

int counter =0;

while(counter<=num_of_iterations)
    {
        std::cout<<"n = "<<counter<<std::endl;
    std::cout<<"i"<<std::setw(15)<<"U"<<std::endl;

    std::cout<<0<<std::setw(15)<<K<<std::endl;
    for (int i = 0; i<m+1; i++)
        {
        std::cout<<i+1<<std::setw(15)<<u[i]<<std::endl;
        }
    u_bar[0] = (1/alpha(0))*(u[0] - gamma(0)*K - beta(0)*u[1]);

    if (g[0] >= (omega*u_bar[0]+(1-omega)*u[0]) )
        {
        u[0] = g[0];
        }
    else
        {
        u[0] = omega*u_bar[0]+(1-omega)*u[0];
        }


    for(int i = 1; i<m-1; i++)
        {
        u_bar[i] = (1.0/alpha(i))*( u[i] - gamma(i)*u[i-1] - beta(i)*u[i+1] );

        if (g[i] >= (omega*u_bar[i]+(1-omega)*u[i]) )
            {
            u[i] = g[i];
            }
        else
            {
            u[i] = omega*u_bar[i]+(1-omega)*u[i];
            }
        }

    u_bar[m-1] = (1/alpha(m-1))*(u[m-1] - gamma(m-1)*u[m-2]);

    if (g[m-1] >= (omega*u_bar[m-1]+(1-omega)*u[m-1]) )
        {
        u[m-1] = g[m-1];
        }
    else
        {
        u[m-1] = omega*u_bar[m-1]+(1-omega)*u[m-1];
        }

    counter ++;

    }

return 0;
}
