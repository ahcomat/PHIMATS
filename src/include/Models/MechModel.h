/**
 * @file MechModel.h
 * @brief A class based on mechanical models to interface with PETSc based global
 * stiffness matrix `A`, solution vector `x` and RHS `a`.
 * 
 * @details Main functions:
 * - Calls `BaseElemMech::CalcElemStiffMatx` to build the local stiffness matrix.
 * - Assembles the global stiffness matrix.
 * - Initializes and applies boundary conditions.
 * - Manages the output by writing to H5file_out.
 * 
 * @todo - Consider inheritance.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-05-24
 * 
 * @copyright Copyright (C) 2024 Abdelrahman Hussein
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

#ifndef MECHMODEL_H
#define MECHMODEL_H

#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include "BaseModel.h"
#include "H5IO.h"
#include "FiniteElements/Mechanics/BaseElemMech.h"
#include "Materials/Mechanics/BaseMechanics.h"

class MechModel: public BaseModel{

public:

MechModel(H5IO& H5File_in);
~MechModel();

/**
 * @brief Set nodal quantities `nodStres`, `nodStran`, `nodCount` and `Fint` to zeros.
 * 
 */
void setZeroNodVals();

/**
 * @brief Initializes and preallocates the RHS `vecFext`, solution `x` and 
 *        stiffness matrix `A.
 * 
 * @param elements 
 */
void InitializePETSc(vector<BaseElemMech*> elements);

/**
 * @brief Set Dirichlet boundary conditions in the RHS `vecFext` and global stiffness matrix `A`.
 * 
 */
void setDirichBC();

/**
 * @brief Updates displacement increment.
 * 
 */
void UpdateDisp();

/**
 * @brief Get the number of steps to apply load.
 * 
 * @return int 
 */
int get_nSteps() const;

/**
 * @brief Calculates the element stiffness matrix. 
 * 
 * @param elements 
 * @param mats 
 */
void CalcElemStiffMatx(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 */
void Assemble(vector<BaseElemMech*> elements);

/**
 * @brief Reads and initializes Dirichlet BCs.
 * 
 * @param H5File_in 
 */
void InitializeDirichBC(H5IO& H5File_in);

/**
 * @brief Pass reference of RHS (to solver).
 * 
 * @return Vec& 
 */
Vec& getB();

/**
 * @brief Pass a reference of the solution vector (to solver).
 * 
 * @return Vec& 
 */
Vec& getX();

/**
 * @brief Pass a reference fo the stiffness matrix (to solver).
 * 
 * @return Mat& 
 */
Mat& getA();

/**
 * @brief Calculates the Fint, strains and stresses. Also Calculates the stress nodal values.
 * 
 * @param elements 
 */
void CalcStres(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out, const string iStep);

private:

/// @brief Tolerance for Newton-Raphson iterations.
const double tol = 1e-6; 

/// @brief Maximum number of iterations.
const int max_iter = 10;  

/// @brief Frequency of updating the stiffness matrix. 
const int NR_freq = 3;      

/// @brief Counter for NR iterations.
int iterCounter = 0;      

/// @brief Number of element displacement dofs.
int nElDispDofs;    

/// @brief Nodal stress.
T_nodStres nodStres;

/// @brief Nodal strain.
T_nodStres nodStran;  

/// @brief Counter for integration points surrounding nodes.
vector<int> nodCount;     

// PETSc ------------------------

/// @brief For calculating the internal force vector.
PetscReal* Fint = NULL;     
/// @brief Indices for `VecSetValues`.
PetscInt* indices = NULL;       
/// @brief Prescribed   
PetscReal* presZeros = NULL;    
/// @brief L2-Norm 
PetscReal l2norm;

/// @brief External force vector.
Vec vecFext; 

/// @brief Residual vector.
Vec vecR; 

/// @brief solution vector.         
Vec vecDisp;

/// @brief Displacement increment vector.         
Vec vecDeltaDisp;

/// @brief solution vector.         
Vec vecFint;

/// @brief The global coefficient (stiffness) matrix.
Mat matA;     

/// @brief `SNES` object.         
SNES snes;  

};
#endif