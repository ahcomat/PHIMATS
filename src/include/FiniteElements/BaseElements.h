/**
 * @file BaseElements.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for elements. BaseElements are designed to behave as
 *        a container for element-set. 
 * @date 2024-05-22
 *  
 * @todo Consider having only elements rather than element set-container.
 *       But this will require another class as element container. 
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

#ifndef BASEELEMENTS_H
#define BASEELEMENTS_H

#include <vector>

#include "Matrix.h"
#include "H5IO.h"

using namespace std;

class BaseElements{

public:

virtual ~BaseElements() = default;

/**
 * @brief Get the total number of DOfs for element set.
 * 
 * @return int 
 */
int get_nDof() const { return nDof; };

/**
 * @brief Get the number of elements. 
 * 
 * @return int 
 */
int get_nElements() const { return nElements; };

/**
 * @brief Returns const reference to the vector of element stiffness matrix `k_ll`.
 * 
 * @return const vector<T_ElStiffMatx>& 
 */
const T_ElStiffMatx& getElStiffMatx() const { return elStiffMatxVariant; }

protected:

/// @brief Total number of elements in element set.
int nElements;         

/// @brief Total number of nodes for element set.
int nNodes;      

/// @brief Total number of DOFs for element set.
int nDof;

/// @brief Logger object for handeling interface messages.
Logger& logger;

/// @brief Int-pt volume.
vector<vector<double>> intPtVol;

/// @brief Gauss points in natural coordinates. 
vector<vector<double>> gaussPts;  

/// @brief Node connectivity.
vector<vector<int>> elemNodeConn;  

/// @brief Global element IDs. 
vector<int> elemIDs;                   

/// @brief Variant for returning elStiffMatx. Passed to `Model` for assembling global stiffness matrix.
T_ElStiffMatx elStiffMatxVariant;      

};
#endif