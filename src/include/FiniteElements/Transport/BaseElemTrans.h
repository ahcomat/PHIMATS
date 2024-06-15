/**
 * @file BaseElemTrans.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for transport elements, i.e. with only concentration/temperature DOFs. 
 * @date 2024-06-13
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASEELEMTRANS_H
#define BASEELEMTRANS_H

#include "FiniteElements/BaseElements.h"

class BaseElemTrans: public BaseElements{

public:

BaseElemTrans(int nElDim, int nElNodes, int nElConDofs, int nElGauss)
    : nElDim(nElDim), nElNodes(nElNodes), nElConDofs(nElConDofs), nElGauss(nElGauss) {};

/**
 * @brief Get dimensions of the element. 
 * 
 * @return int 
 */
int get_nDim() const { return nElDim; };

/**
 * @brief Get the number of element con DOFs. Same as number of nodes. 
 * 
 * @return int 
 */
int get_nElConDofs() const { return nElConDofs; };

/**
 * @brief Return a const reference to the `elemConDof`. 
 * 
 * @return const vector<vector<int>>& 
 */
const vector<vector<int>>& get_elemConDof() const { return elemConDof; };

/**
 * @brief Reads the data `nElements`, `nElementSets` and `elemNodeConn` from hdf5 file.
 * 
 * @param H5File_in  
 */
void ReadElementsData(H5IO &H5File_in, int iSet);

/**
 * @brief Calculates the element stiffness matrix.
 */
virtual void CalcElemStiffMatx(T_DMatx DMatx, double s) = 0;

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
virtual void CalcFlux(T_DMatx KMatx, const double* globalBuffer, T_nodFlux& nodFlux, vector<int>& nodCount) = 0;

protected:

const int nElDim;              /// @brief Spatial dimensions of the element.
const int nElNodes;            /// @brief Number of nodes per element.
const int nElConDofs;          /// @brief Number of element concentration (temperature) dofs.
const int nElGauss;            /// @brief Number of gauss points.

vector<vector<int>> elemConDof;    /// @brief Element concentration (temperature) dofs. In this case, it is identical to `elemNodeConn`.

};
#endif