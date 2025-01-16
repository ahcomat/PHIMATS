/**
 * @file BaseElemMech.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for mechanics elements, i.e. with only displacement DOFs.
 * @date 2024-05-22
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

#ifndef BASEELEMMECH_H
#define BASEELEMMECH_H

#include "FiniteElements/BaseElements.h"
#include "Materials/Mechanics/BaseMechanics.h"

class BaseElemMech: public BaseElements{

public:

BaseElemMech(int nElDim, int nElNodes, int dispDofs, int nElStres, int nElDispDofs, int nElGauss2, string matModel)
    : nElDim(nElDim), nElNodes(nElNodes), dispDofs(dispDofs), nElStres(nElStres),
      nElDispDofs(nElDispDofs), nElGauss(nElGauss2), materialModel(matModel) {};

/**
 * @brief Get dimensions of the element. 
 * 
 * @return int 
 */
int get_nDim() const { return nElDim; };

/**
 * @brief Get the number of element disp DOFs. 
 * 
 * @return int 
 */
int get_nElDispDofs() const { return nElDispDofs; };

/**
 * @brief Return a const reference to the `elemDispDof`. 
 * 
 * @return const vector<vector<int>>& 
 */
const vector<vector<int>>& get_elemDispDof() const { return elemDispDof; };

/**
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in  
 */
void ReadElementsData(H5IO &H5File_in, int iSet);

/**
 * @brief Returns the displacement dofs associated with element `iElem`.
 * 
 * @param iElem 
 * @param dispDof 
 */
void CalcElemDispDof(int iElem, vector<int>& dispDof);

/**
 * @brief Calculates the element stiffness matrix.
 */
virtual void CalcElemStiffMatx(T_DMatx DMatx) = 0;

/**
 * @brief Calculates the `Fint`, strains and stresses. Also evaluates the stress nodal values.
 * 
 * @param DMatx 
 * @param x 
 * @param globalBuffer 
 * @param nodStresFlag 
 * 
 */
virtual void CalcStres(T_DMatx DMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount) = 0;

/**
 * @brief Function for mapping int-point values to nodal points.
 * 
 * @param nodStres 
 * @param nodStran 
 * @param nodStran_e 
 * @param nodStran_p 
 * @param nodStran_eq 
 * @param nodStres_eq 
 * @param nodCount 
 */
virtual void CalcNodVals(T_nodStres& nodStres, T_nodStres& nodStran,T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<int>& nodCount) = 0;

/**
 * @brief Calculates the int-pt strain increment.
 * 
 * @param globalBuffer 
 */
virtual void CalcElDStran(const double* globalBuffer) = 0;

/**
 * @brief Calculates the int-pt total strain.
 * 
 * @param globalBuffer 
 */
virtual void CalcElStran(const double* globalBuffer) = 0;

/**
 * @brief Calls return mapping for plasticity.
 * 
 * @param mat Plastic material object.
 * @param updateStiffMat Flag for updating element stiffness matrix.
 * @param iStep Current step.
 */
virtual void CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep) = 0;

/**
 * @brief Updated _old values.
 * 
 */
virtual void getNew() = 0;

/**
 * @brief Calculates the internal force vector.
 * 
 * @param Fint Buffer for vecFint.
 */
virtual void CalcFint(double* Fint) = 0;

protected:

/// @brief Spatial dimensions of the element.
const int nElDim; 

/// @brief Number of nodes per element.        
const int nElNodes; 

/// @brief Number of displacement dofs per node. 
const int dispDofs;  

/// @brief Stress/strain components.
const int nElStres;      

/// @brief Number of element displacement dofs.
const int nElDispDofs;    

/// @brief Number of gauss points.
const int nElGauss;       

/// @brief Flag for material model [`Elastic`, `ElastoPlastic`]
const string materialModel;     

/// @brief Element displacement dofs.
vector<vector<int>> elemDispDof;    

};
#endif