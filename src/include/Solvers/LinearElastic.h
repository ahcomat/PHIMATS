#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include <petsc.h>

class LinearElastic{

public:

LinearElastic(Mat &A);

~LinearElastic();

/**
 * @brief Solve the linear system `Ax=b`
 * 
 * @param A 
 * @param x 
 * @param b 
 */
void Solve(Vec &x, Vec &b);

private:

KSP ksp;

PC pc;

};
#endif