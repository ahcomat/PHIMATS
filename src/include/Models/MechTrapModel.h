/**
 * @file MechTrapModel.h
 * @brief To be removed.
 * 
 * @details Main functions:
 * - Calls `BaseElemTrap::CalcElemStiffMatx` to build the local diffusivity and capacitance matrix.
 * - Assembles the global capacitance and conductivity matrix.
 * - Initializes and applies boundary conditions.
 * - Manages the output by writing to H5file_out.
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

#ifndef MECHTRAPMODEL_H
#define MECHTRAPMODEL_H

#include <petscvec.h>
#include <petscmat.h>
#include "BaseModel.h"
#include "FiniteElements/Trapping/BaseElemTrap.h"
#include "Materials/Trapping/BaseTrapping.h"
#include "H5IO.h"

class MechTrapModel: public BaseModel{

public:

MechTrapModel(vector<BaseElemTrap*> elements, H5IO& H5File_in);
~MechTrapModel();

/**
 * @brief Set nodal quantities `nodFlux` to zeros.
 * 
 */
void setZero_nodFlux();

/**
 * @brief Initializes and preallocates the RHS `F`, solution `x` and 
 *        conductivity matrix `K_D` and capacitance matrix `C`.
 * 
 * @param elements 
 */
void InitializePETSc(vector<BaseElemTrap*> elements);

/**
 * @brief Calculates the gradients of phi and writes them in `H5File_out` 
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteGradSigmaH(vector<BaseElemTrap*> elements, H5IO& H5File_out);

/**
 * @brief Calculates the element stiffness matrix. 
 * 
 * @param elements 
 * @param mats 
 */
void CalcElemStiffMatx(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 */
void Assemble(vector<BaseElemTrap*> elements);

/**
 * @brief Reads and initializes Dirichlet BCs.
 * 
 * @param H5File_in 
 */
void InitializeDirichBC(H5IO& H5File_in);

void UpdateFT(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats);

/**
 * @brief Set Dirichlet boundary conditions in the RHS `b` and global stiffness matrix `A`.
 * 
 */
void setDirichBC();

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
 * @brief Calculates the flux. Also Calculates the maps the flux to nodal values.
 * 
 * @param elements 
 */
void CalcFlux(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats);

/**
 * @brief Write nodal values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(vector<BaseElemTrap*> elements, H5IO &H5File_out, const string iStep);

private:

int nElConDofs;     /// @brief Number of element concentration (temp) dofs.
double dt;          /// @brief Time increment.
double T;           /// @brief Temperature.

double* FTBuffer;  /// @brief buffer array for PETSc data

T_nodStres nodFlux;          /// @brief Nodal flux.
vector<double> nodCount;     /// @brief Counter for integration points surrounding nodes.

// PETSc ------------------------

Vec F;      /// @brief RHS vector.
Vec x;      /// @brief solution vector.
Mat K_D;    /// @brief The global diffusivity matrix.
Mat MKT;    /// @brief The global [M-KT] matrix.

};
#endif