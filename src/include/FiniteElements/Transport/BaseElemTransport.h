/**
 * @file BaseElemTransport.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for transport elements, i.e. with only concentration/temperature DOFs. 
 * @date 2024-06-13
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

#ifndef BaseElemTransport_H
#define BaseElemTransport_H

#include "FiniteElements/BaseElements.h"

class BaseElemTransport: public BaseElements{

public:

BaseElemTransport(int nElDim, int nElNodes, int nElConDofs, int nElGauss)
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
 * @brief Return const reference to the vector of element capacitance matrix c_ii.
 * 
 * @return const vector<T_ElStiffMatx>& 
 */
const T_ElStiffMatx& getElCapMatx() const { return elCapMatxVariant; }

/**
 * @brief Evaluates the int-pt flux vector. Also evaluates the flux at the nodes.
 * 
 * @param KMatx 
 * @param globalBuffer 
 * @param nodFlux 
 * @param nodCount 
 */
virtual void CalcFlux(T_DMatx KMatx, const double* globalBuffer, T_nodStres& nodFlux, vector<double>& nodCount) = 0;

/**
 * @brief Calculates the volume averaged concentration
 * 
 * @param globalBuffer Solution vector
 * @return double Volume averaged concentration
 */
virtual double CalcAvCon(const double* globalBuffer) = 0;

protected:

/// @brief Number of element concentration (temperature) dofs.
const int nElConDofs;          

/// @brief Time increment.
double dt;                          

/// @brief Element concentration (temperature) dofs. In this case, it is identical to `elemNodeConn`.
vector<vector<int>> elemConDof;    

/// @brief Variant for returning elCapMatx. 
T_ElStiffMatx elCapMatxVariant;   

};
#endif