/**
 * @file LinearTransport.cxx
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief A wrapper for PETSc linear solver `ksp` mainly for linear heat and mass transport.
 * @date 2024-06-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LINEARTRANSPORT_H
#define LINEARTRANSPORT_H

#include "BaseSolver.h"

class LinearTransport: public BaseSolver{

public:

LinearTransport(Mat &A);

~LinearTransport() override;

void UpdateKSP(Mat &A);

/**
 * @brief Solve the linear system `Ax=b`.
 * 
 * @param A 
 * @param x 
 * @param b 
 */
void SolveTransport(Vec &x, Vec &F);

};
#endif