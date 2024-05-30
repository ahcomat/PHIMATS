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
 * Updates (when, what and who)
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
 * @brief Get the total number of DOfs.
 * 
 * @return int 
 */
int get_nTotDof() const { return nTotDof; };

/**
 * @brief Get the number of elements. 
 * 
 * @return int 
 */
int get_nElements() const { return nElements; };

/**
 * @brief Evaluates the element stiffness matrix.
 */
virtual void CalcElemStiffMatx(T_DMatx DMatx) = 0;

/**
 * @brief Return const reference to the vector of element stiffness matrix k_ll.
 * 
 * @return const vector<T_ElStiffMatx>& 
 */
const T_ElStiffMatx& getElStiffMatx() const { return elStiffMatxVariant; }

/**
 * @brief Write element specific int-pts output averaged over the nodes.
 * 
 * @param H5File_out 
 */
virtual void WriteOut(H5IO &H5File_out) = 0;

protected:

int nElementSets;   /// @brief Number of element sets
int nElements;      /// @brief Total number of elements.
int nNodes;         /// @brief Total number of nodes.
int nTotDof;        /// @brief Total number of DOFs.
int nPresDofs;      /// @brief Number of prescribed displacement dofs.

vector<vector<double>> gaussPts;    /// @brief Gauss points in natural coordinates. 
vector<vector<int>> elemNodeConn;   /// @brief Node connectivity.

T_ElStiffMatx elStiffMatxVariant;   /// @brief Variant for returning elStiffMatx. 

};
#endif