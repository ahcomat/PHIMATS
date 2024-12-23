/**
 * @file LinearElastic.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear elasticity.
 * @date 2024-05-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LINEARELASTIC_H
#define LINEARELASTIC_H

#include "BaseSolver.h"

class LinearElastic: public BaseSolver{

public:

LinearElastic(Mat &A);

~LinearElastic() override;

/**
 * @brief Solve the linear system `Ax=b`.
 * 
 * @param A 
 * @param x 
 * @param b 
 */
void Solve(Vec &x, Vec &b);

};
#endif