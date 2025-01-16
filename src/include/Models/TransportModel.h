/**
 * @file TransportModel.h
 * @brief A class based on heat and mass transport models to interface with PETSc based global
 * capacitance matrix `C`, conductivity matrix `K`, solution vector `x` and RHS `F`.
 * 
 * @details Main functions:
 * - Calls `BaseElemTransport::CalcElemStiffMatx` to build the local diffusivity and capacitance matrix.
 * - Assembles the global capacitance and conductivity matrix.
 * - Initializes and applies boundary conditions.
 * - Manages the output by writing to H5file_out.
 * 
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @date 2024-06-15
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

#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include <petscvec.h>
#include <petscmat.h>
#include "BaseModel.h"
#include "H5IO.h"
#include "FiniteElements/Transport/BaseElemTransport.h"
#include "Materials/Transport/BaseTransport.h"

class TransportModel: public BaseModel{

public:

TransportModel(vector<BaseElemTransport*> elements, H5IO& H5File_in);
~TransportModel();

/**
 * @brief Set nodal quantities `nodStres`, `nodStran` and `nodCount` to zeros.
 * 
 */
void setZero_nodFlux();

/**
 * @brief Initializes and preallocates the RHS `F`, solution `x` and 
 *        conductivity matrix `K` and capacitance matrix `C`.
 * 
 * @param elements 
 */
void InitializePETSc(vector<BaseElemTransport*> elements);

/**
 * @brief Calculates the element stiffness matrix. 
 * 
 * @param elements 
 * @param mats 
 */
void CalcElemStiffMatx(vector<BaseElemTransport*> elements, vector<BaseTransport*> mats);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 */
void Assemble(vector<BaseElemTransport*> elements);

/**
 * @brief Reads and initializes Dirichlet BCs.
 * 
 * @param H5File_in 
 */
void InitializeDirichBC(H5IO& H5File_in);

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
void CalcFlux(vector<BaseElemTransport*> elements, vector<BaseTransport*> mats);

/**
 * @brief Write the total concentration of the time step
 * 
 * @param H5File_out 
 * @param tStep 
 */
void WriteAvCon(vector<BaseElemTransport*> elements, H5IO &H5File_out, const int iStep);

/**
 * @brief Writes the average exit flux in the x direction.
 * 
 * @param H5File_out 
 * @param tStep 
 */
void WriteAvFlux(H5IO &H5File_out, const int tStep);

/**
 * @brief Write nodal values of concentration.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(H5IO &H5File_out, const string iStep);

private:

int nElConDofs;     /// @brief Number of element concentration (temp) dofs.
double dt;          /// @brief Time increment.
int nExitNodes;     /// @brief Number of exit nodes.

T_nodStres nodFlux;          /// @brief Nodal flux.
vector<double> nodCount;     /// @brief Counter for integration points surrounding nodes.

vector<int> ExitNodeIDs;     /// @brief Exit nodes IDs.

// PETSc ------------------------

Vec F;   /// @brief RHS vector.
Vec x;   /// @brief solution vector.
Mat K;   /// @brief The global conductivity matrix.
Mat M;   /// @brief The global capacitance matrix.

};
#endif