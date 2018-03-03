#include <iostream>
#include <cmath>

// Function for fx
double fx(double x)
{
    // Define variable f
    double f;

    // Function f(x) from equation (3)
    f = tanh(100.0*(x-2.0))*0.1*(3.0-x);

    f = f*sqrt(x*x - x*sqrt(2.0) + 0.5005);

    // Returns f
    return f;
}

// Function for GQ1
double GQ1(double a, double b)
{
    // Define variables
    double w1 = 2.0, integral;

    // Apply rule to arbitray integration interval
    integral = 0.5*(b-a)*w1*fx(0.5*(a+b));

    // Returns the integral
    return integral;
}

// Function for GQ3
double GQ3(double a, double b)
{
    // Define variables
    int i;
    double wi[] = {5.0/9, 8.0/9, 5.0/9};
    double x0 = sqrt(3.0/5);
    double xi[3] = {-x0, 0, x0};
    double integral = 0.0, dx;

    // Define dx
    dx = 0.5*(b-a);

    for (i=0;i<3;i++)
    {
        // Apply rule to arbitray integration interval
        integral = integral + wi[i]*fx(xi[i]*dx + 0.5*(a+b));
    }

    integral = integral*dx;

    // Returns integral
    return integral;
}

double AdaptInt(double a, double b, double tol)
{
    // Define variables
    int iter, N = 1, N1, n=1000;
    double *xl, *dx, *err, *tol1;
    double xk, I1, I3, sum1, sum3, err_tot;

    xl= new double[n];
    dx= new double[n];
    err=new double[n];
    tol1=new double[n];

    //Initialize the first interval
    xl[0] = a;
    xl[1] = b;
    dx[0] = b-a;

    err_tot = 1e3;

    // While total estimate of error > tolerance
    while (err_tot > tol)
    {
        iter++;

        err_tot = 0.0;
        I1 = 0.0;
        I3 = 0.0;

        for (int i = 0; i < N; i++)
        {
            // t_k=((x_k-x_(k-1))/b-a))*t
            tol1[i] = tol*dx[i]/(b-a);

            xk = xl[i] + dx[i];

            // Total of GQ1^(k)
            sum1 = GQ1(xl[i], xk);

            // Total of GQ3^(k)
            sum3 = GQ3(xl[i], xk);

            // Estimate of error
            err[i] = fabs(sum3-sum1);
         // Total of estimate of error
            err_tot = err_tot + err[i];

            I1 = I1 + sum1;
            I3 = I3 + sum3;
        }

        // Output iteration, N, GQ1, GQ3 and error
        std::cout<<"  "<< N << " " << I1 <<" " << I3 << " " << err_tot <<std::endl;

        N1 = N;

        for (int i = 0; i < N1; i++)
        {
            // If estimate of error > t_k
            if (err[i] > tol1[i])
            {
                N1++;

                // Shift all xl and dx from, i+1 to next index
                for (int j = N1; j > i+1; j--)
                {
                    xl[j] = xl[j-1];
                }

                for (int j = N1-1; j > i+1; j--)
                {
                    dx[j] = dx[j-1];
                    err[j] = err[j-1];
                    tol1[j] = tol1[j-1];
                }

                // Add new point
                xl[i+1] = 0.5*(xl[i] + xl[i+1]);
                dx[i] = 0.5*dx[i];
                dx[i+1] = dx[i];
                i++;
            }
        }

        N = N1;
    }

    delete [] xl;
    delete [] err;
    delete [] dx;
    delete [] tol1;

    return I1;
}

// Main function
int main(int argc, char* argv[])
{
    // Define variables
    int N;
    double a = 0.0, b = 3.0, sum1, sum2;
    double tol = 1.0e-6;

    sum1 = AdaptInt(a,b,tol);
    std::cout<<"sum 1 = "<<sum1<<"\n";
    return sum1;
}
