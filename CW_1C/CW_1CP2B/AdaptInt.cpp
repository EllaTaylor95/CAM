#include <iostream>
#include <cmath>
#include <iomanip>

double CGQ1(int n, double a, double b);
double CGQ3(int n, double a, double b);

double adapt_int(double a, double b, double tau)
{
int n = 1; //number of  intervals
int n_new;
int N = 1000; //max number of mesh points
int counter = 0;
double Est= 10;
double Gauss_quad1 = 0.0;
double Gauss_quad3 = 0.0;


double* tau_k;
double* Est_k;
double* x;

x = new double[N];
tau_k = new double[N];
Est_k = new double[N];

//Setting end points for single interval
x[0] = a;
x[1] = b;

std::cout<< " i "<<std::setw(10)<<" N "<<std::setw(20)<<" CGQ1 "<<std::setw(30)<<" CGQ3 "<<std::setw(40)<<" Sum Est_k "<<std::endl;

while (Est > tau)
    {
    counter = counter + 1;
    //Est =/= 0 for while loop to begin. Now recalculating Est
    //Resetting Gauss_Quad calculations to zero
    Est = 0.0;
    Gauss_quad1 = 0.0;
    Gauss_quad3 = 0.0;

    for (int k = 0; k<n; k++)
        {
        std::cout<<"x["<<k<<"] = "<<x[k]<<"\n";
        //Storing each tolerance
        tau_k[k] = tau*(x[k+1] - x[k])/(b-a);
        //Storing each error
        Est_k[k] = fabs(CGQ3(1,x[k],x[k+1]) - CGQ1(1,x[k],x[k+1]));

        Est = Est + Est_k[k];

        Gauss_quad1 = Gauss_quad1 + CGQ1(1, x[k],x[k+1]);
        Gauss_quad3 = Gauss_quad3 + CGQ3(1, x[k],x[k+1]);
        }

    std::cout << counter  <<std::setw(10)<<n<<std::setw(20)<< Gauss_quad1 <<std::setw(30)<<Gauss_quad3 <<std::setw(40)<<Est<<std::endl;

    n_new = n;

    for(int k = 0; k<n_new; k++)
        {
        if(Est_k[k] > tau_k[k])
            {
            n_new= n_new+1;

            //Shifting entries of x along to make space for new mesh point
            for(int i=n_new; i>k+1; i--)
                {
                x[i] = x[i-1];
                }
            for(int i=n_new-1; i>k+1 ; i--)
                {
                tau_k[i] = tau_k[i-1];
                Est_k[i] = Est_k[i-1];
                }
            //New mesh point
            x[k+1] = 0.5*(x[k]+x[k+1]);
            k=k+1;
            }
        }

    n=n_new;
    }

delete[] x;
delete[] tau_k;
delete[] Est_k;

return Gauss_quad1;
}
