#ifndef "CW2_P4HEADERDEF"
#define "CW2_P4HEADERDEF"

#include <cmath>

class ParabolicODE

public:
    //Constructor
    ParabolicODE(const int m);

    //Destructor
    ~ParabolicODE();

    //Accessor
    double* GetRHS() const;

    //Calulating error
    void CalculateErr(double* u, const double* u_actual) const;

    //Returns solution for tridiagonal matrix eq
    void TridiagonalMatrixSolver() const;

    double normal_CFD(double x);


private:

    //Cannot remember what goes on in this section here

    //Default constructor
    ParabolicODE();

    double* mpRHS;

    double mStepSize;

    double mDeltaT;


#endif
