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

BaseElemMech(int nElDim, int nElNodes, int nElGauss, int dispDofs, int nElStres, int nElDispDofs,  string matModel, Logger& logger)
    : BaseElements(nElDim, nElNodes, nElGauss, logger), dispDofs(dispDofs), 
      nElStres(nElStres), nElDispDofs(nElDispDofs), materialModel(matModel) {};

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
virtual void CalcElemStiffMatx(T_DMatx CMatx) = 0;

 /**
  * @brief Calculates strains, stresses and internal force for linear elastic material. Maps strain/stress to nodal values.
  * 
  * @param CMatx Stiffness matrix
  * @param globalBuffer Displacement solution vector
  * @param Fint Internal force vector
  * @param nodStres Stresse mapped to the nodes
  * @param nodStran Strains mapped to the nodes
  * @param nodCount 
  */
virtual void CalcStres(T_DMatx CMatx, const double* globalBuffer, double* Fint, T_nodStres& nodStres, T_nodStres& nodStran, vector<int>& nodCount) = 0;

/**
 * @brief Maps int-point values to nodes for output (plasticity).
 * 
 * @param nodStres Stress
 * @param nodStran Total strain
 * @param nodStran_e Elastic strain
 * @param nodStran_p Plastic strain
 * @param nodStran_eq Equivalent plastic strain
 * @param nodStres_eq Equivalent (Von Mises) stress
 * @param nodStres_h  Hydrostatic stress
 * @param nodCount 
 */
virtual void CalcNodVals(T_nodStres& nodStres, T_nodStres& nodStran,T_nodStres& nodStran_e, T_nodStres& nodStran_p, vector<double>& nodStran_eq, vector<double>& nodStres_eq, vector<double>& nodStres_h, vector<double>& nodRho, vector<int>& nodCount) = 0;

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

/// @brief Number of displacement dofs per node. 
const int dispDofs; 

/// @brief Stress/strain components.
const int nElStres;                 

/// @brief Number of element displacement dofs.
const int nElDispDofs;

/// @brief Flag for material model [`Elastic`, `ElastoPlastic`]
const string materialModel;     

/// @brief Element displacement dofs.
vector<vector<int>> elemDispDof;    

/// @brief Int-pt equivalent plastic strain [nElStres].
vector<vector<double>> elStran_eq;

/// @brief Int-pt equivalent plastic strain (last convergged increment) [nElStres].
vector<vector<double>> elStran_eq_old;
    
/// @brief Int-pt equivalent stress (von Mises) [nElStres].
vector<vector<double>> elStres_eq;

/// @brief Int-pt hydrostatic stress [nElStres].
vector<vector<double>> elStres_h;

/// @brief Int-pt normalized dislocation density [nElStres].
vector<vector<double>> elRho;

};
#endif