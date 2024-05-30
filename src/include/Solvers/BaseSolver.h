/**
 * @file BaseSolver.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear elasticity.
 * @date 2024-05-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASESOLVER_H
#define BASESOLVER_H

#include <petscvec.h>

class BaseSolver{

public:

virtual ~BaseSolver() = default;

/**
 * @brief Solve the linear system `Ax=b`.
 * 
 * @param A 
 * @param x 
 * @param b 
 */
virtual void Solve(Vec &x, Vec &b) = 0;

};
#endif