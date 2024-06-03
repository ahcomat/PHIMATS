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

#include <petscsys.h>
#include <petscvec.h>

#include "FiniteElements/BaseElements.h"

using namespace std;

class BaseElemMech: public BaseElements{

public:

BaseElemMech(int nDim, int nElNodes, int dispDofs, int nStres, int nElDispDofs, int nGauss)
    : nDim(nDim), nElNodes(nElNodes), dispDofs(dispDofs), nStres(nStres),
      nElDispDofs(nElDispDofs), nGauss(nGauss) {};

/**
 * @brief Get dimensions of the element. 
 * 
 * @return int 
 */
int get_nDim() const { return nDim; };

/**
 * @brief Get the number of element disp DOFs. 
 * 
 * @return int 
 */
int get_nElDispDofs() const { return nElDispDofs; };

/**
 * @brief Return a const reference to the disp DOFs of element `iElem`. 
 * 
 * @return const vector<vector<int>>& 
 */
const vector<int>& get_elemDispDof(int iElem) const { return elemDispDof.at(iElem); };

/**
 * @brief Calculates the `Fint`, strains and stresses. Also evaluates the stress nodal values 
 *        if `nodStresFlag=true`.
 * 
 * @param DMatx 
 * @param x 
 * @param globalBuffer 
 * @param nodStresFlag 
 * 
 * @todo 
 * - Remove `nodStresFlag`.
 */
virtual void CalcStres(T_DMatx DMatx, const double* globalBuffer) = 0;

protected:

const int nDim;           /// @brief Spatial dimensions of the element.
const int nElNodes;       /// @brief Number of nodes per element.
const int dispDofs;       /// @brief Number of displacement dofs. 
const int nStres;         /// @brief Stress/strain components.
const int nElDispDofs;    /// @brief Number of element displacement dofs.
const int nGauss;         /// @brief Number of gauss points.

vector<vector<int>> elemDispDof;    /// @brief Element displacement dofs.

};
#endif