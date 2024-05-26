/**
 * @file Elements.h
 * @author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
 * @brief The base class for elements. Elements are designed to behave as
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

#ifndef ELEMENTS_H
#define ELEMENTS_H

#include<vector>

#include"Matrix.h"

using namespace std;

class Elements{

public:
Elements(int nDim, int nElNodes, int dispDofs, int nStres, int nElDispDofs, int nGauss)
    : nDim(nDim), nElNodes(nElNodes), dispDofs(dispDofs), nStres(nStres),
      nElDispDofs(nElDispDofs), nGauss(nGauss) {};

virtual ~Elements() = default;

/**
 * @brief Evaluates the element stiffness matrix.
 */
virtual void CalcElemStiffMatx(T_DMatx DMatx) = 0;

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
 * @brief Return a const reference to the `elStiffMatx`.
 * 
 * @return const vector<T_elStiffMatx>& 
 */
virtual const T_elStiffMatx& get_elStiffMatx() const = 0;

protected:

const int nDim;           /// @brief Spatial dimensions of the element.
const int nElNodes;       /// @brief Number of nodes per element.
const int dispDofs;       /// @brief Number of displacement dofs. 
const int nStres;         /// @brief Stress/strain components.
const int nElDispDofs;    /// @brief Number of element displacement dofs.
const int nGauss;         /// @brief Number of gauss points.

int nElementSets;   /// @brief Number of element sets
int nElements;      /// @brief Total number of elements.
int nNodes;         /// @brief Total number of nodes.
int nTotDof;        /// @brief Total number of DOFs.
int nPresDofs;      /// @brief Number of prescribed displacement dofs.

vector<vector<double>> gaussPts;    /// @brief Gauss points in natural coordinates. 
vector<vector<int>> elemNodeConn;   /// @brief Node connectivity.
vector<vector<int>> elemDispDof;    /// @brief Element displacement dofs.

T_elStiffMatx elStiffMatxVariant;   /// @brief Variant for returning elStiffMatx. 

};
#endif