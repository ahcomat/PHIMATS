/**
 * @file BaseElements.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for elements. BaseElements are designed to behave as
 *        a container for element-set. 
 * @date 2024-05-22
 * 
 * @copyright Copyright (c) 2024
 * 
 * @todo Consider having only elements rather than element set-container.
 *       But this will require another class as element container. 
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
 * @brief Calculates the element stiffness matrix.
 */
virtual void CalcElemStiffMatx(T_DMatx DMatx) = 0;

/**
 * @brief Return const reference to the vector of element stiffness matrix k_ll.
 * 
 * @return const vector<T_ElStiffMatx>& 
 */
const T_ElStiffMatx& getElStiffMatx() const { return elStiffMatxVariant; }

protected:

int nElements;         /// @brief Total number of elements in element set.
int nNodes;            /// @brief Total number of nodes for element set.
int nDof;              /// @brief Total number of DOFs for element set.
// int nPresDofs;      /// @brief Number of prescribed displacement dofs.


vector<vector<double>> gaussPts;    /// @brief Gauss points in natural coordinates. 
vector<vector<int>> elemNodeConn;   /// @brief Node connectivity.
vector<int> elemIDs;                   /// @brief Global element IDs. 

T_ElStiffMatx elStiffMatxVariant;   /// @brief Variant for returning elStiffMatx. 

};
#endif