/**
 * @file BaseModel.h
 * @brief A base class for `Models`.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-06-20
 * 
 * @copyright Copyright (C) 2025 Abdelrahman Hussein
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
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
PetscInt nPresDofs = 0;            /// @brief number of prescribed dofs.
PetscInt* presDofs = NULL;     /// @brief Array to hold the prescribed dofs.
PetscScalar* presVals = NULL;  /// @brief Array to hold the prescribed values.

};
#endif