/**
 * @file TrappingModel.h
 * @brief A class based on heat and mass transport models to interface with PETSc based global
 * capacitance matrix `C`, conductivity matrix `K_D`, solution vector `x` and RHS `F`.
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

#ifndef TrappingModel_H
#define TrappingModel_H

#include <petscvec.h>
#include <petscmat.h>
#include "BaseModel.h"
#include "FiniteElements/Trapping/BaseElemTrap.h"
#include "Materials/Trapping/BaseTrapping.h"
#include "H5IO.h"

class TrappingModel: public BaseModel{

public:

TrappingModel(vector<BaseElemTrap*> elements, H5IO& H5File_in);
~TrappingModel();

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
void WriteGradPhi(vector<BaseElemTrap*> elements, H5IO& H5File_out);

/**
 * @brief Updates the temperature according the the heating rate HR.
 * 
 * @param iStep Time step
 * @param HR Heating rate
 */
void UpdateTemp(const int iStep, double HR);

/**
 * @brief Updates time increment `dt`.
 * 
 * @param elements Element sets vector
 * @param dtNew New time increment
 */
void Update_dt(vector<BaseElemTrap*> elements, double dtNew);

/**
 * @brief Write the current temperature.
 * 
 * @param H5File_out Output hdf5 file
 * @param tStep Time step
 */
void WriteTemp(H5IO &H5File_out, const int iStep);

/**
 * @brief Calculates the element stiffness matrix.
 * 
 * @param elements 
 * @param mats 
 * @param updateTemp 
 */
void CalcElemStiffMatx(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats);

/**
 * @brief Assemble the global stiffness matrix.
 * 
 * @param elements 
 * @param updateTemp 
 */
void Assemble(vector<BaseElemTrap*> elements, bool updateTemp=false);

/**
 * @brief Reads initial conditions from H5File.
 * 
 * @param H5File 
 * @param iStep 
 */
void ReadInitialCon(H5IO& H5File, const int iStep);

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
 * @brief Function for updating `F`.
 * 
 */
void Update_F();

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

void WriteInPtCoords(vector<BaseElemTrap*> elements, H5IO &H5File_out);

/**
 * @brief Calculates the flux. Also Calculates the maps the flux to nodal values.
 * 
 * @param elements 
 */
void CalcFlux(vector<BaseElemTrap*> elements, vector<BaseTrapping*> mats);

/**
 * @brief Write nodal flux field values.
 * 
 * @param H5File_out 
 * @param iStep 
 */
void WriteFlux(H5IO &H5File_out, const string iStep);

void WriteIntPtFlux(H5IO &H5File_out, const string iStep);

/**
 * @brief Write nodal concentration values.
 * 
 * @param elements 
 * @param H5File_out 
 */
void WriteOut(H5IO &H5File_out, const string iStep);

/**
 * @brief Write the volume averaged concentration of the time step.
 * 
 * @param H5File_out Ouput hdf5 file handle
 * @param tStep Time step
 */
void WriteAvCon(vector<BaseElemTrap*> elements, H5IO &H5File_out, const int iStep);

/**
 * @brief Writes the average exit flux in the x direction.
 * 
 * @param H5File_out 
 * @param tStep 
 */
void WriteExitFlux(H5IO &H5File_out, const int tStep);

private:

/// @brief Number of element concentration (temp) dofs.
int nElConDofs;     

/// @brief Number of exit nodes.
int nExitNodes;

/// @brief Time increment.
double dt;    

/// @brief Initial temperature.
double T0;

/// @brief Current temperature.
double T;

/// @brief Total time.
double TotTime = 0;     

/// @brief Total number of integration points.
int  nTotGuasPts = 0;   

/// @brief int-pt flux.
T_nodStres intPtFlux;  

/// @brief Nodal flux.
T_nodStres nodFlux;    

/// @brief Counter for integration points surrounding nodes.
vector<double> nodCount;     

/// @brief Exit nodes IDs.
vector<int> ExitNodeIDs;     

// PETSc ------------------------

/// @brief RHS vector.
Vec vecF;  

/// @brief solution vector.
Vec vecx;    

/// @brief The global stiffness matrix.
Mat matK;      

/// @brief The global [M-KT] matrix.
Mat matM;

};
#endif