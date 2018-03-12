#include <iostream>
#include <cmath>
#include <iomanip>

void funct(double* f, double* y, double alpha, double beta)
{
f[0] = y[1];
f[1] = -pow(alpha,2)*sin(y[0])-pow(beta,2)*y[1];
}

void Newton_Method(double TOL, double* y_1, double* y_2, int n_max, double alpha, double beta, double h)
{
double* u;
u = new double[2];

int counter = 0;

double *f;
f = new double[2];
double F[2];

double F_norm = 10.0;
double determinant;
double Energy;

std::cout<< "n"<<std::setw(7)<<"k"<<std::setw(12)<<"y1"<<std::setw(18)<<"y2"<<std::setw(25)<<"F_norm"<<std::endl;

for (int i = 1; i <= n_max; i++)
    {
    //Setting initial guess y_n+1(0) = y_n

   std::cout<<std::endl<<i;//<<std::endl;

    y_1[i] = y_1[i-1];
    y_2[i] = y_2[i-1];

    // calculating initial norm
     u[0] = y_1[i];
     u[1] = y_2[i];
     funct(f,u,alpha,beta);

    F[0] = y_1[i] - y_1[i-1] - h*f[0];
    F[1] = y_2[i] - y_2[i-1] - h*f[1];
    F_norm = sqrt(pow(F[0],2) + pow(F[1],2));
    //std::cout<<"F[1] = "<<F[1]<<"\n";

    while (F_norm >= TOL)
        {
        F_norm = 0;

     //   std::cout<<std::setw(7)<<counter;
        counter++;

        // update y1 and y2
        u[0] = y_1[i];
        u[1] = y_2[i];
        funct(f,u,alpha,beta);
        determinant = (1.0 + pow(h,2)*pow(alpha,2)*cos(y_1[i]));

        y_1[i] +=  (-1.0/determinant) * (F[0] + h*F[1]);
    //   std::cout<<std::setw(12)<<y_1[i];
        y_2[i] += (-1.0/determinant) * ( (-h*pow(alpha,2)*cos(y_1[i])*F[0]) + F[1]);
    //    std::cout<<std::setw(18)<<y_2[i];

        // calculate norm for updated y and y2 values
        u[0] = y_1[i];
        u[1] = y_2[i];
        funct(f,u,alpha,beta);
        F[0] = y_1[i] - y_1[i-1] - h*f[0];
        //std::cout<<"F[0] = "<<F[0]<<"\n";
        F[1] = y_2[i] - y_2[i-1] - h*f[1];
        //std::cout<<"F[1] = "<<F[1]<<"\n";
        F_norm = sqrt(pow(F[0],2) + pow(F[1],2));
        //y_1[i]<<"  " << y_2[i] << "\n";

      //  std::cout<<std::setw(25)<<F_norm<<std::endl;
        }
        counter=0;

        Energy = 0.5*pow(y_2[i],2)+pow(alpha,2)*(1-cos(y_1[i]));
        std::cout<<std::setw(25)<<Energy<<std::endl;

    }
delete[] f;
delete[] u;
}

int main()
{
double alpha = 2.0;
double beta = 0.0;
double pi = 3.14159265359;

double T = 8.0;
int n_max = 32;
double h = T/(double)(n_max);

double* y_1;
double* y_2;

y_1 = new double[n_max+1];
y_2 = new double[n_max+1];


y_1[0] = pi/2.0;
y_2[0] = 0.0;

for (int i = 1; i<n_max; i++)
{
    y_1[i] = 0;
    y_2[i] = 0;
}

double TOL = pow(10,-12);

//From backward Euler have F(y_n+1) = y_n+1 - y_n -hf(y_n+1,t+1) = 0
//Then apply Newtons method to this

Newton_Method(TOL, y_1, y_2, n_max, alpha, beta, h);



delete[] y_1;
delete[] y_2;

return 0;
}
