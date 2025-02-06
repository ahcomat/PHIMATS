/**
 * @file MechModel.h
 * @brief A class for mechanical models. Interface with `PETSc`.
 * 
 * @details Main functions:
 * - Calls `BaseElemMech::CalcElemStiffMatx` to build the local stiffness matrix.
 * - Assembles the global stiffness matrix.
 * - Initializes and applies boundary conditions.
 * - Manages non-linear solvers through `SNES`.
 * - Manages the output by writing to H5file_out.
 * 
 * @todo - Consider inheritance.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-05-24
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

/**
 * @brief Constructor. Reads main data from hdf5 file. 
 * 
 * @param H5File_in Input hdf5 file. 
 * @param NR_update Frequency of updating the Jacobian and its predonditioner.
 */

/**
 * @brief Construct a new Mech Model object
 * 
 * @param H5File_in 
 */
MechModel(H5IO& H5File_in, const int NR_update = 3);
~MechModel();

/**
 * @brief Reads and initializes Dirichlet BCs.
 * 
 * @param H5File_in Input hdf5 file. 
 */
void InitializeDirichBC(H5IO& H5File_in);

/**
 * @brief Set Dirichlet boundary conditions in the RHS `vecFext` and global stiffness matrix `A`.
 * 
 */
void setDirichBC();

/**
 * @brief Initializes and preallocates PETSc objects.
 * 
 * @param elements Vector of elements.
 */
void InitializePETSc(vector<BaseElemMech*> elements);

/**
 * @brief Set nodal quantities `nodStres`, `nodStran`, `nodCount` and `Fint` to zeros.
 * 
 */
void setZeroNodVals();

/**
 * @brief Get the number of steps to apply load increment.
 * 
 * @return int Number of steps.
 */
int get_nSteps() const;

/**
 * @brief Calculates the element stiffness matrix. 
 * 
 * @param elements Vector of elements.
 * @param mats Vector of materials.
 */
void CalcElemStiffMatx(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats);

/**
 * @brief Helper function for `Assemble`.
 * 
 * @param elStiffMatx_ptr 
 * @param elemDispDof_ptr 
 * @param nElDispDofs 
 * @param nElements 
 * @return PetscErrorCode 
 */
PetscErrorCode AssembleElementMatrix(const auto* elStiffMatx_ptr,
                                     const std::vector<std::vector<int>>& elemDispDof_ptr,
                                     PetscInt nElDispDofs,
                                     PetscInt nElements);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 * @param elements Vector of elements.
 * @return PetscErrorCode.
 */
PetscErrorCode Assemble(vector<BaseElemMech*> elements);

/**
 * @brief SNES solve of the current step.
 * 
 * @param elements Elements vector.
 * @param mats Materials vector. 
 * @param iStep Current step. 
 */
void SolveSNES(vector<BaseElemMech*> elements, vector<BaseMechanics*> mats, int iStep);

/**
 * @brief Casts `CalcResidual` to static to pass the residual to the SNES solver. 
 * 
 * @param snes SNES solver.
 * @param deltaU Solution vector.
 * @param R Residual vector
 * @param ctx Application context.
 * @return PetscErrorCode 
 */
static PetscErrorCode ResidualCallback(SNES snes, Vec deltaU, Vec R, void *ctx);

/**
 * @brief Casts `Assemble` to static to pass the Jacobian to the SNES solver.
 * 
 * @param snes SNES solver.
 * @param deltaU Solution vector.
 * @param J Jacobian matix.
 * @param P Preconditoiner matrix. Same as J.
 * @param ctx Application context.
 * @return PetscErrorCode 
 */
static PetscErrorCode JacobianCallback(SNES snes, Vec deltaU, Mat J, Mat P, void *ctx);

/**
 * @brief Calculates the residual vector `vecR`.
 * 
 * @param deltaU Displacement increment (solution vector).
 * @param elements Vector of elements.
 * @param mats Materials vector.
 * @param iStep Current step.
 * @return PetscErrorCode 
 */
PetscErrorCode CalcResidual(Vec deltaU, Vec R, vector<BaseElemMech*> elements, vector<BaseMechanics*> mats, int iStep);

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
 * @brief Maps int-point values to nodes for output.
 * 
 * @param elements 
 */
void CalcNodVals(vector<BaseElemMech*> elements);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(vector<BaseElemMech*> elements, H5IO &H5File_out, const string iStep);

private:

/// @brief Maximum number of iterations.
const int max_iter = 10;  

/// @brief Frequency of updating the stiffness matrix. 
const int NR_freq;      

/// @brief Counter for NR iterations.
int iterCounter = 0;  

///@brief Flag for updating stiffness matrix
bool updateStiffMat = true;

/// @brief Number of element displacement dofs.
int nElDispDofs;    

/// @brief Nodal stress.
T_nodStres nodStres;

/// @brief Nodal strain.
T_nodStres nodStran;  

/// @brief Nodal elastic strain.
T_nodStres nodStran_e;

/// @brief Nodal plastic strain.
T_nodStres nodStran_p;

/// @brief Nodal equivalent plastic strain.
vector<double> nodStran_eq;

/// @brief Nodal equivalent stress (Von Mises).
vector<double> nodStres_eq;

/// @brief Nodal hydrostatic stress.
vector<double> nodStres_h;

/// @brief Counter for integration points surrounding nodes.
vector<int> nodCount;     

// PETSc ------------------------

/// @brief Buffer for calculating the internal force vector.
PetscReal* Fint = NULL;     

/// @brief Indices for `VecSetValues`.
PetscInt* indices = NULL;   

/// @brief To set the Residual prescribed dofs to zeros.    
PetscReal* presZeros = NULL;    

/// @brief External force vector.
Vec vecFext; 

/// @brief Residual vector.
Vec vecR; 

/// @brief Displacement vector.         
Vec vecDisp;

/// @brief Displacement increment (solution) vector.         
Vec vecDeltaDisp;

/// @brief Internal force vector.         
Vec vecFint;

/// @brief The global coefficient (stiffness) matrix.
Mat matA;

/// @brief `KSP` object for controlling the linear system.
KSP ksp;

/// @brief `PC` object for controlling the preconditioner. 
PC pc;

/// @brief `SNES` object.         
SNES snes; 

/// @brief `SNES` convergence reason.
SNESConvergedReason reason; // Convergence reason

/// @brief Error code.
PetscErrorCode ierr;        // Error code

/// @brief Define application context structure for SNES. 
struct AppCtx {
    vector<BaseElemMech*> elements;  // Elements vector
    vector<BaseMechanics*> mats;     // Material vector
    int iStep;                       // Current step 
    MechModel *mechModel;            // Pointer to <this>
};

};
#endif