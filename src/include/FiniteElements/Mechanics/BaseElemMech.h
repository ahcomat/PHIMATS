/**
 * @file BaseElemMech.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for mechanics elements, i.e. with only displacement DOFs.
 * @date 2024-05-22
 * 
 * @copyright Copyright (c) 2024 
 *  
 */

#ifndef BASEELEMMECH_H
#define BASEELEMMECH_H

#include "FiniteElements/BaseElements.h"
#include "Materials/Mechanics/BaseMechanics.h"

class BaseElemMech: public BaseElements{

public:

BaseElemMech(int nElDim, int nElNodes, int dispDofs, int nElStres, int nElDispDofs, int nElGauss2)
    : nElDim(nElDim), nElNodes(nElNodes), dispDofs(dispDofs), nElStres(nElStres),
      nElDispDofs(nElDispDofs), nElGauss(nElGauss2) {};

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
 * @brief Calculates the int-pt total strain.
 * 
 * @param globalBuffer 
 */
virtual void CalcElStran(const double* globalBuffer) = 0;

/**
 * @brief Calls return mapping for plasticity.
 * 
 * @param mats 
 */
virtual void CalcRetrunMapping(BaseMechanics* mat, const bool& updateStiffMat, int iStep) = 0;

virtual void CalcFint(double* Fint) = 0;

protected:

const int nElDim;         /// @brief Spatial dimensions of the element.
const int nElNodes;       /// @brief Number of nodes per element.
const int dispDofs;       /// @brief Number of displacement dofs per node. 
const int nElStres;       /// @brief Stress/strain components.
const int nElDispDofs;    /// @brief Number of element displacement dofs.
const int nElGauss;       /// @brief Number of gauss points.

string materialModel;     /// @brief Flag for material model [`Elastic`, `ElasoPlastic`]

vector<vector<int>> elemDispDof;    /// @brief Element displacement dofs.

};
#endif