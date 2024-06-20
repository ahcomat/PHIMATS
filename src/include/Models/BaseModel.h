/**
 * @file BaseModel.h
 * @brief A base class for `Models`.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-06-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASEODEL_H
#define BASEODEL_H

#include <petscsys.h>

class BaseModel{

protected:

int nElementSets;   /// @brief Number of element sets
int nTotNodes;      /// @brief Total number of nodes.
int nTotDofs;       /// @brief Total number of DOFs.
int nTotElements;   /// @brief Total number of elements.
int nDim;           /// @brief Spatial dimensions of the model.
int nElements;      /// @brief Number of elements per element set.
int nSteps;         /// @brief Number of steps to apply the load.

// PETSc ------------------------

const double* globalBuffer;  /// @brief buffer array for PETSc data

// Boundary conditions
PetscInt nPresDofs;            /// @brief number of prescribed dofs.
PetscInt* presDofs = NULL;     /// @brief Array to hold the prescribed dofs.
PetscScalar* presVals = NULL;  /// @brief Array to hold the prescribed values.

};
#endif