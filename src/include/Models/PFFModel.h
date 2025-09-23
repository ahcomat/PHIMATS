/**
 * @file PFFModel.h
 * @brief A class handel phase-field fracture models. Interfaces with `PETSc`.
 * 
 * @details Main functions:
 * - Evaluates field dependent work energy densities. 
 * - Calculates PFF driving forces.
 * - Calls `BaseElemPFF::CalcElemStiffMatx`.
 * - Assembles the global stiffness matrix.
 * - Manages the output by writing to H5file_out.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2025-09-21
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

#ifndef PFFMode_H
#define PFFMode_H

#include <petscvec.h>
#include <petscmat.h>
#include "BaseModel.h"
#include "Logger.h"
#include "FiniteElements/PFF/BaseElemPFF.h"
#include "FiniteElements/Mechanics/BaseElemMech.h"
#include "Materials/PFF/BasePFF.h"
#include "Materials/Mechanics/BaseMechanics.h"
#include "H5IO.h"

class PFFModel: public BaseModel{

public:

PFFModel(vector<BaseElemPFF*> elements, H5IO& H5File_in, Logger& logger);
~PFFModel();

/**
 * @brief Initializes and preallocates the RHS `F`, solution `x` and 
 *        conductivity matrix `K_D` and capacitance matrix `C`.
 * 
 * @param elements 
 */
void InitializePETSc(vector<BaseElemPFF*> pffElem);

/**
 * @brief Sets `pffElem->elem_wc` to constant value `wc`.
 * 
 * @param pffElem PFF elements vector
 * @param pffMat PFF materials vector
 */
void set_const_wc(vector<BaseElemPFF*> pffElem, vector<BasePFF*> pffMat);

/**
 * @brief Calculates the spectral decomposition of the strain energy density. 
 * 
 * @param pffElem PFF elements vector
 * @param mechElem Mechanics elements vector
 * @param mechMat Mechanics materials vector
 */
void CalcPsiSpectral(vector<BaseElemPFF*> pffElem, vector<BaseElemMech*> mechElem, vector<BaseMechanics*> mechMat);

/**
 * @brief Calculates elastic crack driving force based on `psi_plus`.
 * 
 * @param pffElem PFF elements vector
 */
void CalcDrivFrocElas(vector<BaseElemPFF*> pffElem);

/**
 * @brief Calculates the element stiffness matrix.
 * 
 * @param pffElem PFF elements vector
 */
void CalcElemStiffMatx(vector<BaseElemPFF*> pffElem);

/**
 * @brief Helper function for `Assemble`.
 * 
 * @param elMatx_ptr 
 * @param elemConDof_ptr 
 * @param nElConDofs 
 * @param nElements 
 * @param globalMat 
 */
void AssembleElementMatrix(const auto* elMatx_ptr,
                           const vector<vector<int>>& elemPhiDof_ptr,
                           PetscInt nElPhiDofs,
                           PetscInt nElements,
                           Mat globalMat);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 * @param pffElem Vector of element set pointers.
 * @param assembleM Flag for updating the capacitance matrix, default=true.
 */
void Assemble(vector<BaseElemPFF*> pffElem);

void CalcFp(vector<BaseElemPFF*> pffElem);

void Calc_gPhi_d(vector<BaseElemPFF*> pffElem);

/**
 * @brief Get the number of steps to apply load.
 * 
 * @return int 
 */
int get_nSteps() const;

/**
 * @brief Pass reference of RHS (to solver).
 * 
 * @return Vec& 
 */
Vec& getF();

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
Mat& getK();

/**
 * @brief Write nodal concentration values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(H5IO &H5File_out, const string iStep);

private:

/**
 * @brief Logger object for handeling interface messages.
 * 
 */
Logger& logger;

/// @brief Number of element phi dofs.
int nElPhiDofs = 0;     

/// @brief Total number of integration points.
int  nTotGuasPts = 0;        

// PETSc ------------------------

/// @brief Buffer for calculating the RHS vector.
PetscReal* Fp = NULL;  

/// @brief Indices for `VecSetValues`.
PetscInt* indices = NULL;   

/// @brief RHS vector.
Vec vecFp;  

/// @brief solution vector.
Vec vecx;    

/// @brief solution vector.
Vec vecPhi;  

/// @brief The global stiffness matrix.
Mat matK;      

};
#endif