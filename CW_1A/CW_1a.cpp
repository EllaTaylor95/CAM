#include <iostream>

//Writing a code that implements the tridiagonal matrix algorithm for solving an nxn matrix

void tridiagonal_matrix_solver(int n, double* x, double* lower, double* diag, double* upper, double* f)
    {
    //Elimination stage

    //f[0] = f[0] and d[0] = d[0] no change
    for (int i = 1; i<n; i++)
        {
        diag[i] = diag[i] - ((upper[i-1]*lower[i-1])/diag[i-1]);
        f[i]    = f[i] - ((f[i-1]*lower[i-1])/diag[i-1]);
        }

    //Backsolving
    //Bottom row is a special case
    x[n-1] = f[n-1]/diag[n-1];

    for (int i = n-2; i >= 0; i--)
        {
        x[i] = (f[i] - upper[i]*x[i+1])/diag[i];
        }

    }
