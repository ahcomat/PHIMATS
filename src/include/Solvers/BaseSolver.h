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

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
class BaseSolver{

public:

virtual ~BaseSolver() = default;

protected:

KSP ksp;        /// @brief `KSP` object.
PC pc;          /// @brief Pre-conditioner.

};
#endif