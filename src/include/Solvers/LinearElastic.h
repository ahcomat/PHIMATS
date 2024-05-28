/**
 * @file LinearElastic.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear elasticity.
 * @date 2024-05-28
 * 
 * @copyright Copyright (c) 2024
 * 
 * Updates (when, what and who)
 * 
 */

#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include <petscsys.h>
#include <petscmat.h>
#include <petscksp.h>

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

/**
 * @brief `KSP` object.
 * 
 */
KSP ksp;

/**
 * @brief Pre-conditioner.
 * 
 */
PC pc;

};
#endif